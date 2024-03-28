import os
import tempfile
import shutil
import multiprocessing
import pickle
from copy import deepcopy
from enum import Enum
from typing import Optional, List, Tuple, Dict, Any
from typing_extensions import Literal

import rdkit.Chem as Chem

from ccdc.docking import Docker as DockerGold
from ccdc.io import MoleculeReader, EntryWriter
from ccdc.protein import Protein
from pydantic import BaseModel

from dockstream.core.Schrodinger.Glide_docker import Parallelization
from dockstream.utils.enums.logging_enums import LoggingConfigEnum
from dockstream.utils.execute_external.Gold import GoldExecutor

from dockstream.core.docker import Docker
from dockstream.core.Gold.Gold_result_parser import GoldResultParser
from dockstream.utils.enums.Gold_enums import GoldLigandPreparationEnum
from dockstream.utils.enums.Gold_enums import GoldTargetKeywordEnum, GoldExecutablesEnum, GoldOutputEnum
from dockstream.utils.general_utils import gen_temp_file

from dockstream.utils.translations.molecule_translator import MoleculeTranslator
from dockstream.utils.dockstream_exceptions import DockingRunFailed


class GoldFitnessFunction(str, Enum):
    GOLDSCORE = "goldscore"
    CHEMSCORE = "chemscore"
    ASP = "asp"
    PLP = "plp"


class GoldResponseValue(str, Enum):
    FITNESS = "fitness"
    VALUE = "value"


class GoldParameters(BaseModel):
    prefix_execution: Optional[str] = None
    binary_location: Optional[str] = None
    receptor_paths: Optional[List[str]] = None
    time_limit_per_compound: Optional[int] = None
    parallelization: Optional[Parallelization]
    fitness_function: GoldFitnessFunction
    response_value: GoldResponseValue = GoldResponseValue.FITNESS
    early_termination: bool
    autoscale: float  # Autoscale percentage. very fast: 10, medium: 50, very slow: 100.
    ndocks: int = 10
    diverse_solutions: Optional[Tuple[bool, Optional[int], Optional[float]]] = None   # If diverse solutions is enabled this will be (True, cluster size, rmsd), otherwise (False, None, None). TODO: rework for GUI.

    def get(self, key: str) -> Any:
        """Temporary method to support nested_get"""
        return self.dict()[key]

_LP = GoldLigandPreparationEnum()
_TK = GoldTargetKeywordEnum()
_EE = GoldExecutablesEnum()
_ROE = GoldOutputEnum()
_LE = LoggingConfigEnum()


class Gold(Docker):
    """Interface to the Gold backend."""

    backend: Literal["Gold"] = "Gold"
    parameters: GoldParameters

    _target_dict: Dict = None
    _Gold_executor: GoldExecutor = None
    _scoring_function_parameters: Dict[str, str] = None

    class Config:
        underscore_attrs_are_private = True

    def __init__(self, **run_parameters):
        # invoke base class's constructor first
        super().__init__(**run_parameters)

        # prepare and check Gold backend availability
        self._check_Gold_backend_availability()

        # parse the fitness function and response value set
        self._parse_fitness_function()

        # set the tag name for the scoring function and whether minimial or maximum values are better
        self._scoring_function_parameters = self._get_scoring_function_parameters()

    def _check_Gold_backend_availability(self):
        self._Gold_executor = GoldExecutor(
            prefix_execution=self.parameters.prefix_execution,
            binary_location=self.parameters.binary_location)
        if not self._Gold_executor.is_available():
            raise DockingRunFailed("Cannot initialize Gold docker, as Gold backend is not available - abort.")
        self._logger.log(f"Checked Gold backend availability (prefix_execution={self.parameters.prefix_execution}).", _LE.DEBUG)

    def _parse_fitness_function(self):
        self._logger.log(f"Set fitness function to {self.parameters.fitness_function} and response value to {self.parameters.response_value}.", _LE.DEBUG)

    def _initialize_cavity(self, settings):
        # load the target dictionary specification and initialize the cavity
        target_path = self.parameters.receptor_paths[0]
        with open(target_path, "rb") as file:
            self._target_dict = pickle.load(file)
            self._logger.log(f"Loaded pickled cavity dictionary stored in file {target_path}.", _LE.DEBUG)
            if self._target_dict[_TK.VERSION] != _TK.CURRENT_VERSION:
                self._logger.log(f"Version of pickled target ({self._target_dict[_TK.VERSION]}) is not the same as DockStream's ({_TK.CURRENT_VERSION}).", _LE.WARNING)
        self._logger.log(f"Unpacked the target dictionary.", _LE.DEBUG)

        tmpdir = tempfile.mkdtemp()
        if self._target_dict[_TK.CAVITY_METHOD] == _TK.CAVITY_METHOD_REFERENCE:
            # write ligand to temporary file (ending copied over in settings)
            tmp_ref_ligand_path = gen_temp_file(suffix=self._target_dict[_TK.REFERENCE_LIGAND_FILENAME], dir=tmpdir)
            with open(tmp_ref_ligand_path, 'w') as file:
                for line in self._target_dict[_TK.REFERENCE_LIGAND]:
                    file.write(line)
                self._logger.log(f"Wrote temporary ligand file {tmp_ref_ligand_path} with {len(self._target_dict[_TK.REFERENCE_LIGAND])} lines.", _LE.DEBUG)

            # write target PDB to temporary file
            tmp_target_path = gen_temp_file(suffix=".pdb", dir=tmpdir)
            with open(tmp_target_path, 'w') as file:
                for line in self._target_dict[_TK.TARGET_PDB]:
                    file.write(line)
                self._logger.log(f"Wrote temporary target file {tmp_target_path} with {len(self._target_dict[_TK.TARGET_PDB])} lines.", _LE.DEBUG)

            # build the cavity
            ref_ligand = MoleculeReader(filename=tmp_ref_ligand_path)[0]
            self._prepare_protein(settings, tmp_target_path)
            protein = settings.proteins[0]
            settings.binding_site = settings.BindingSiteFromLigand(protein,
                                                                   ref_ligand,
                                                                   distance=self._target_dict[_TK.CAVITY_REFERENCE_DISTANCE])
            settings.reference_ligand_file = tmp_ref_ligand_path
        elif self._target_dict[_TK.CAVITY_METHOD] == _TK.CAVITY_METHOD_POINT:
            raise NotImplementedError
            # origin (x,x,x)
            # distance x
        else:
            raise DockingRunFailed("Specified cavity determination method not defined for GOLD.")
        self._logger.log(f"Initialized GOLD Protein.BindingSite with method {self._target_dict[_TK.CAVITY_METHOD]}.", _LE.DEBUG)

    def add_molecules(self, molecules: list):
        """This method overrides the parent class, docker.py add_molecules method. This method appends prepared
        ligands to a list for subsequent docking. Note, that while internally we will store the ligands for "GOLD"
        in RDkit format, they will need to be written out as an SDF file before docking can commence later.

        :param molecules: A list that is to contain all prepared ligands for subsequent docking
        :type molecules: list
        :raises NotImplementedError: Each backend must override the parent class, docker.py add_molecules method.
            Inability to do so or a bug causing incorrect implementation will raise a NotImplementedError
        """
        mol_trans = MoleculeTranslator(self.ligands, force_mol_type=_LP.TYPE_RDKIT)
        mol_trans.add_molecules(molecules)
        self.ligands = mol_trans.get_as_rdkit()
        self._docking_performed = False

    def _generate_temporary_input_output_files(self, start_indices, sublists):
        # in case singletons are handed over, wrap them in a list for "zipping" later
        if not isinstance(start_indices, list):
            start_indices = [start_indices]
        if not isinstance(sublists, list):
            sublists = [sublists]

        tmp_output_dirs = []
        tmp_input_sdf_paths = []
        tmp_output_sdf_paths = []
        for start_index, sublist in zip(start_indices, sublists):
            # generate temporary input files and output directory
            cur_tmp_output_dir = tempfile.mkdtemp()
            cur_tmp_sdf = gen_temp_file(prefix=str(start_index), suffix=".sdf", dir=cur_tmp_output_dir)
            
            # write-out the temporary input file
            writer = Chem.SDWriter(cur_tmp_sdf)
            one_written = False
            for ligand in sublist:
                # initialize all ligands (as they could have failed)
                if ligand.get_molecule() is not None:
                    mol = deepcopy(ligand.get_molecule())
                    mol.SetProp("_Name", ligand.get_identifier())
                    one_written = True
                    writer.write(mol)
            writer.close()
            if one_written is False:
                if os.path.isdir(cur_tmp_output_dir):
                    shutil.rmtree(cur_tmp_output_dir)
                continue

            # add the path to which "_dock_subjob()" will write the result SDF
            output_sdf_path = gen_temp_file(prefix=str(start_index), suffix="_result.sdf", dir=cur_tmp_output_dir)
            tmp_output_dirs.append(cur_tmp_output_dir)
            tmp_output_sdf_paths.append(output_sdf_path)
            tmp_input_sdf_paths.append(cur_tmp_sdf)
        return tmp_output_dirs, tmp_input_sdf_paths, tmp_output_sdf_paths

    def _dock(self, number_cores: int):
        # partition ligands into sublists and distribute to processor cores for docking
        start_indices, sublists = self.get_sublists_for_docking(number_cores=number_cores)
        number_sublists = len(sublists)
        self._logger.log(f"Split ligands into {number_sublists} sublists for docking.", _LE.DEBUG)
        sublists_submitted = 0
        slices_per_iteration = min(number_cores, number_sublists)

        while sublists_submitted < len(sublists):
            upper_bound_slice = min((sublists_submitted + slices_per_iteration), len(sublists))
            cur_slice_start_indices = start_indices[sublists_submitted:upper_bound_slice]
            cur_slice_sublists = sublists[sublists_submitted:upper_bound_slice]

            # generate paths and initialize molecules (so that if they fail, this can be covered)
            tmp_output_dirs, tmp_input_sdf_paths, \
            tmp_output_sdf_paths = self._generate_temporary_input_output_files(cur_slice_start_indices,
                                                                               cur_slice_sublists)

            # run in parallel; wait for all subjobs to finish before proceeding
            processes = []
            for chunk_index in range(len(tmp_output_dirs)):
                p = multiprocessing.Process(target=self._dock_subjob, args=(tmp_input_sdf_paths[chunk_index],
                                                                            tmp_output_sdf_paths[chunk_index],
                                                                            tmp_output_dirs[chunk_index]))
                processes.append(p)
                p.start()
            for p in processes:
                p.join()

            # add the number of input sublists rather than the output temporary folders to account for cases where
            # entire sublists failed to produce an input structure
            sublists_submitted += len(cur_slice_sublists)

            # load the chunks and recombine the result; add conformations
            for chunk_index in range(len(tmp_output_dirs)):
                # this is a protection against the case where empty (file size == 0 bytes) files are generated due to
                # a failure during docking
                if not os.path.isfile(tmp_output_sdf_paths[chunk_index]) or os.path.getsize(tmp_output_sdf_paths[chunk_index]) == 0:
                    continue

                for molecule in Chem.SDMolSupplier(tmp_output_sdf_paths[chunk_index], removeHs=False):

                    # it can happen, that ligands have "impossible chemistry" and will be loaded by RDkit as "None"
                    if molecule is None:
                        continue

                    # parse the molecule name (sorted by FITNESS not the score) which looks like:
                    # "0:0|0xa6enezm|sdf|1|dock6"
                    cur_conformer_name = str(molecule.GetProp("_Name")).split(sep='|')[0]

                    # add molecule to the appropriate ligand
                    for ligand in self.ligands:
                        if ligand.get_identifier() == cur_conformer_name:
                            ligand.add_conformer(molecule)
                            break

            # clean-up
            for path in tmp_output_dirs:
                shutil.rmtree(path)
            self._log_docking_progress(number_done=sublists_submitted, number_total=number_sublists)

        # update conformer names to contain the conformer id
        # -> <ligand_number>:<enumeration>:<conformer_number>
        reverse = True if self._get_scoring_function_parameters()[_ROE.BEST] == "max" else False
        for ligand in self.ligands:
            ligand.set_conformers(sorted(ligand.get_conformers(),
                                         key=lambda x: self._get_score_from_conformer(conformer=x),
                                         reverse=reverse))
            ligand.add_tags_to_conformers()

        # log any docking fails
        self._docking_fail_check()

        # generate docking results as dataframe
        result_parser = GoldResultParser(ligands=[ligand.get_clone() for ligand in self.ligands],
                                         fitness_function=self.parameters.fitness_function,
                                         response_value=self.parameters.response_value)
        self._df_results = result_parser.as_dataframe()

        # set docking flag
        self._docking_performed = True

    def _dock_subjob(self, sdf_ligand_path, path_sdf_results, tmp_output_dir):
        # 1) prepare Gold docker: (i) "clone" the docker instance, (ii) set remaining, ligang-specific settings and
        #                         (iii) initialize this chunk's ligands
        cur_docker = DockerGold()
        settings = cur_docker.settings
        settings.output_directory = tmp_output_dir
        settings.output_file = os.path.basename(path_sdf_results)
        settings.output_format = "sdf"
        settings.fitness_function = self.parameters.fitness_function
        settings.early_termination = self.parameters.early_termination
        settings.autoscale = self.parameters.autoscale

        if self.parameters.diverse_solutions is not None:
            settings.diverse_solutions = self.parameters.diverse_solutions

        self._initialize_cavity(settings)

        settings.add_ligand_file(sdf_ligand_path, ndocks=self.parameters.ndocks)

        # 2) write settings file
        settings_file_path = os.path.join(tmp_output_dir, _EE.GOLD_AUTO_CONFIG_NAME)
        settings.write(settings_file_path)
        with open(settings_file_path, 'r') as file:
            self._logger.log(f"Contents of configurations file {settings_file_path}:", _LE.DEBUG)
            for line in file:
                self._logger_blank.log(line.rstrip("\n"), _LE.DEBUG)

        # 3) run Gold docker
        arguments = [settings_file_path]
        execution_result = self._Gold_executor.execute(command=_EE.GOLD_AUTO,
                                                       arguments=arguments,
                                                       check=False)
        self._delay4file_system(path=path_sdf_results)
        self._logger.log(f"Finished sublist (input: {sdf_ligand_path}, output directory: {tmp_output_dir}), with return code '{execution_result.returncode}'.", _LE.DEBUG)

    def _prepare_protein(self, settings, tmp_protein_path):
        protein = Protein.from_file(tmp_protein_path)
        protein.remove_all_waters()
        protein.remove_unknown_atoms()
        protein.add_hydrogens()

        ligands = protein.ligands
        for l in ligands:
            protein.remove_ligand(l.identifier)
        protein_file_name = os.path.join(settings.output_directory, 'clean_%s.mol2' % protein.identifier)

        with EntryWriter(protein_file_name) as writer:
            writer.write(protein)
        settings.add_protein_file(protein_file_name)
        return ligands

    def write_docked_ligands(self, path, mode="all"):
        self._write_docked_ligands(path, mode, mol_type=_LP.TYPE_RDKIT)

    def _get_scoring_function_parameters(self):
        # get the appropriate name of the tag and whether minimal or maximal values are best for
        # the specified scoring function
        if self.parameters.response_value == GoldResponseValue.FITNESS:
            scoring_function_parameters = _ROE.DICT_FITNESS[self.parameters.fitness_function]
        elif self.parameters.response_value == GoldResponseValue.VALUE:
            scoring_function_parameters = _ROE.DICT_VALUE[self.parameters.fitness_function]
        else:
            raise ValueError("Parameter response value must be either fitness or value.")
        self._logger.log(f"Set scoring_function_parameters to {scoring_function_parameters} for obtaining the scores.",
                         _LE.DEBUG)
        return scoring_function_parameters

    def get_scores(self, best_only):
        """This method overrides the parent class, docker.py get_scores method. This method returns a list containing
        all docking scores. This method allows returning the best docking scores only. "best" can mean the minimum
        or maximum values for this given scoring function. By default, it will return the minimum values. Returning
        all the docking scores (of different poses for instance) is also possible if best only is not enforced

        :param best_only: Determines whether the best (either minimum or maximum) docking scores are returned
        :type best_only: boolean, True of False
        :return: list of returned docking scores
        :raises ValueError: If best_only is True but neither "min" nor "max" was specified, a ValueError is raised
        """
        return self._get_scores(best_only=best_only, best=self._scoring_function_parameters[_ROE.BEST])

    def write_result(self, path, mode="all"):
        """This method overrides the parent class, docker.py write_result method.
        This method writes the docking results to a csv file. There is the option to write out either the best
        predicted binding pose per enumeration or all the predicted binding poses. Output for the best predicted
        binding pose per ligand has yet to be implemented

        :param path: Contains information on results output path
        :type path: string
        :param mode: Determines whether the output contains the best predicted binding pose per ligand, the best
            predicted binding pose per enumeration, or all the predicted binding poses
        :type mode: string, optional, default value is "all". Other possible value is "best_per_enumeration"
        :param best: Determines whether lower or higher values are better (typically lower ones)
        :type best: string, optional, default value is "min". Other possible value is "max"
        """
        return self._write_result(path=path, mode=mode, best=self._scoring_function_parameters[_ROE.BEST])

    def _get_score_from_conformer(self, conformer):
        return float(conformer.GetProp(self._scoring_function_parameters[_ROE.TAG]))

    def _sort_conformers(self, conformers: list, best=None) -> list:
        return super()._sort_conformers(conformers=conformers,
                                        best=self._scoring_function_parameters[_ROE.BEST])
