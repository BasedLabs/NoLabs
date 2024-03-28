import os
import tempfile
import shutil
import multiprocessing
from copy import deepcopy
from typing import Optional, List, Any

import rdkit.Chem as Chem
from pydantic import BaseModel, Field
from typing_extensions import Literal

from dockstream.core.Schrodinger.Glide_docker import Parallelization
from dockstream.core.docker import Docker
from dockstream.core.AutodockVina.AutodockVina_result_parser import AutodockResultParser
from dockstream.utils.enums.logging_enums import LoggingConfigEnum
from dockstream.utils.execute_external.AutodockVina import AutodockVinaExecutor
from dockstream.utils.enums.AutodockVina_enums import AutodockVinaExecutablesEnum, AutodockVinaOutputEnum, AutodockResultKeywordsEnum
from dockstream.utils.execute_external.OpenBabel import OpenBabelExecutor
from dockstream.utils.enums.OpenBabel_enums import OpenBabelExecutablesEnum
from dockstream.utils.enums.RDkit_enums import RDkitLigandPreparationEnum
from dockstream.utils.general_utils import gen_temp_file

from dockstream.utils.translations.molecule_translator import MoleculeTranslator
from dockstream.utils.dockstream_exceptions import DockingRunFailed

_LP = RDkitLigandPreparationEnum()
_RKA = AutodockResultKeywordsEnum()
_BEE = OpenBabelExecutablesEnum()
_LE = LoggingConfigEnum()
_ROE = AutodockVinaOutputEnum()
_EE = AutodockVinaExecutablesEnum()


class SearchSpace(BaseModel):
    center_x: float = Field(alias="--center_x")
    center_y: float = Field(alias="--center_y")
    center_z: float = Field(alias="--center_z")
    size_x: float = Field(alias="--size_x")
    size_y: float = Field(alias="--size_y")
    size_z: float = Field(alias="--size_z")

    class Config:
        allow_population_by_field_name = True


class AutodockVinaParameters(BaseModel):
    prefix_execution: Optional[str] = None
    binary_location: Optional[str] = None
    parallelization: Optional[Parallelization]
    receptor_pdbqt_path: Optional[List[str]] = None
    search_space: SearchSpace
    seed: int = 42
    number_poses: int = 1

    def get(self, key: str) -> Any:
        """Temporary method to support nested_get"""
        return self.dict()[key]


class AutodockVina(Docker, BaseModel):
    """Interface to the "AutoDock Vina" backend."""

    backend: Literal["AutoDockVina"] = "AutoDockVina"
    parameters: AutodockVinaParameters

    _ADV_executor: AutodockVinaExecutor = None
    _OpenBabel_executor: OpenBabelExecutor = None

    class Config:
        underscore_attrs_are_private = True

    def __init__(self, **data):
        super().__init__(**data)

    def _initialize_executors(self):
        """Initialize executors and check if they are available."""

        self._ADV_executor = AutodockVinaExecutor(
            prefix_execution=self.parameters.prefix_execution, 
            binary_location=self.parameters.binary_location
        )
        if not self._ADV_executor.is_available():
            raise DockingRunFailed("Cannot initialize AutoDock Vina docker, as backend is not available - abort.")
        self._logger.log(f"Checked AutoDock Vina backend availability (prefix_execution={self.parameters.prefix_execution}).",
                         _LE.DEBUG)

        # initialize the executor for all "OpenBabel" related calls and also check if it is available
        # note, that while there is an "OpenBabel" API (python binding) which we also use, the match to the binary
        # options is not trivial; thus, use command-line here
        self._OpenBabel_executor = OpenBabelExecutor()
        if not self._OpenBabel_executor.is_available():
            raise DockingRunFailed(
                "Cannot initialize OpenBabel external library, which should be part of the environment - abort.")

    def _get_score_from_conformer(self, conformer):
        return float(conformer.GetProp(_RKA.SDF_TAG_SCORE))

    def add_molecules(self, molecules: list):
        """This method overrides the parent class, docker.py add_molecules method. This method appends prepared
        ligands to a list for subsequent docking. Note, that while internally we will store the ligands for
        "AutoDock Vina" in RDkit format, they will need to be written out in AutoDock's pdbqt format for
        subsequent docking

        :param molecules: A list that is to contain all prepared ligands for subsequent docking
        :type molecules: list
        :raises NotImplementedError: Each backend must override the parent class, docker.py add_molecules method.
            Inability to do so or a bug causing incorrect implementation will raise a NotImplementedError
        """
        mol_trans = MoleculeTranslator(self.ligands, force_mol_type=_LP.TYPE_RDKIT)
        mol_trans.add_molecules(molecules)
        self.ligands = mol_trans.get_as_rdkit()
        self._docking_performed = False

    def _write_molecule_to_pdbqt(self, path, molecule) -> bool:
        # generate temporary copy as PDB
        temp_pdb = gen_temp_file(suffix=".pdb")
        Chem.MolToPDBFile(mol=molecule, filename=temp_pdb)

        # Note: In contrast to the target preparation,
        # we will use a tree-based flexibility treatment here -
        # thus, the option "-xr" is NOT used.
        arguments = [temp_pdb,
                     _BEE.OBABEL_OUTPUT_FORMAT_PDBQT,
                     "".join([_BEE.OBABEL_O, path]),
                     _BEE.OBABEL_PARTIALCHARGE, _BEE.OBABEL_PARTIALCHARGE_GASTEIGER]
        self._OpenBabel_executor.execute(command=_BEE.OBABEL,
                                         arguments=arguments,
                                         check=False)

        if os.path.exists(path):
            return True
        else:
            return False

    def _generate_temporary_input_output_files(self, start_indices, sublists):
        # in case singletons are handed over, wrap them in a list for "zipping" later
        if not isinstance(start_indices, list):
            start_indices = [start_indices]
        if not isinstance(sublists, list):
            sublists = [sublists]

        tmp_output_dirs = []
        tmp_input_paths = []
        tmp_output_paths = []
        ligand_identifiers = []
        for start_index, sublist in zip(start_indices, sublists):
            # for "AutoDock Vina", only single molecules can be handled so every sublist is guaranteed at this stage
            # to have only one element
            # TODO: in the course of the redesign of the docking class-backends, also "beautify" this part
            for ligand in sublist:
                # generate temporary input files and output directory
                cur_tmp_output_dir = tempfile.mkdtemp()
                cur_tmp_input_pdbqt = gen_temp_file(prefix=str(start_index), suffix=".pdbqt", dir=cur_tmp_output_dir)
                cur_tmp_output_sdf = gen_temp_file(prefix=str(start_index), suffix=".sdf", dir=cur_tmp_output_dir)

                # write-out the temporary input file
                if ligand.get_molecule() is None:
                    if os.path.isdir(cur_tmp_output_dir):
                        shutil.rmtree(cur_tmp_output_dir)
                    continue
                mol = deepcopy(ligand.get_molecule())
                self._write_molecule_to_pdbqt(cur_tmp_input_pdbqt, mol)
                tmp_output_dirs.append(cur_tmp_output_dir)
                tmp_input_paths.append(cur_tmp_input_pdbqt)
                tmp_output_paths.append(cur_tmp_output_sdf)
                ligand_identifiers.append(ligand.get_identifier())

        return tmp_output_dirs, tmp_input_paths, tmp_output_paths, ligand_identifiers

    def _dock(self, number_cores):

        self._initialize_executors()

        start_indices, sublists = self.get_sublists_for_docking(number_cores=number_cores, enforce_singletons=True)
        number_sublists = len(sublists)
        self._logger.log(f"Split ligands into {number_sublists} sublists for docking.", _LE.DEBUG)

        if not os.path.exists(self.parameters.receptor_pdbqt_path[0]):
            raise DockingRunFailed("Specified PDBQT path to target (receptor) does not exist - abort.")

        sublists_submitted = 0
        slices_per_iteration = min(number_cores, number_sublists)
        while sublists_submitted < len(sublists):
            upper_bound_slice = min((sublists_submitted + slices_per_iteration), len(sublists))
            cur_slice_start_indices = start_indices[sublists_submitted:upper_bound_slice]
            cur_slice_sublists = sublists[sublists_submitted:upper_bound_slice]

            # generate paths and initialize molecules (so that if they fail, this can be covered)
            tmp_output_dirs, tmp_input_paths, tmp_output_paths, \
            ligand_identifiers = self._generate_temporary_input_output_files(cur_slice_start_indices,
                                                                             cur_slice_sublists)

            # run in parallel; wait for all subjobs to finish before proceeding
            processes = []
            for chunk_index in range(len(tmp_output_dirs)):
                p = multiprocessing.Process(target=self._dock_subjob, args=(tmp_input_paths[chunk_index],
                                                                            tmp_output_paths[chunk_index]))
                processes.append(p)
                p.start()
            for p in processes:
                p.join()

            # add the number of input sublists rather than the output temporary folders to account for cases where
            # entire sublists failed to produce an input structure
            sublists_submitted += len(cur_slice_sublists)

            # parse the resulting sdf files
            for path_sdf_results, cur_identifier in zip(tmp_output_paths, ligand_identifiers):
                # add conformations
                if not os.path.isfile(path_sdf_results) or os.path.getsize(path_sdf_results) == 0:
                    continue

                for molecule in Chem.SDMolSupplier(path_sdf_results, removeHs=False):
                    if molecule is None:
                        continue

                    # extract the score from the AutoDock Vina output and update some tags
                    score = self._extract_score_from_VinaResult(molecule=molecule)
                    molecule.SetProp("_Name", cur_identifier)
                    molecule.SetProp(_RKA.SDF_TAG_SCORE, score)
                    molecule.ClearProp(_ROE.REMARK_TAG)

                    # add molecule to the appropriate ligand
                    for ligand in self.ligands:
                        if ligand.get_identifier() == cur_identifier:
                            ligand.add_conformer(molecule)
                            break

            # clean-up
            for path in tmp_output_dirs:
                shutil.rmtree(path)
            self._log_docking_progress(number_done=sublists_submitted, number_total=number_sublists)

        # the conformers are already sorted, but some tags are missing
        # -> <ligand_number>:<enumeration>:<conformer_number>
        for ligand in self.ligands:
            ligand.add_tags_to_conformers()

        # log any docking fails
        self._docking_fail_check()

        # generate docking results as dataframe
        result_parser = AutodockResultParser(ligands=[ligand.get_clone() for ligand in self.ligands])
        self._df_results = result_parser.as_dataframe()

        # set docking flag
        self._docking_performed = True

    def _extract_score_from_VinaResult(self, molecule) -> str:
        result_tag_lines = molecule.GetProp(_ROE.REMARK_TAG).split("\n")
        result_line = [line for line in result_tag_lines if _ROE.RESULT_LINE_IDENTIFIER in line][0]
        parts = result_line.split()
        return parts[_ROE.RESULT_LINE_POS_SCORE]

    def _dock_subjob(self, input_path_pdbqt, output_path_sdf):

        # set up arguments list and execute
        # TODO: support "ensemble docking" - currently, only the first entry is used
        tmp_pdbqt_docked = gen_temp_file(suffix=".pdbqt", dir=os.path.dirname(input_path_pdbqt))
        search_space = self.parameters.search_space
        arguments = [_EE.VINA_RECEPTOR, self.parameters.receptor_pdbqt_path[0],
                     _EE.VINA_LIGAND, input_path_pdbqt,
                     _EE.VINA_CPU, str(1),
                     _EE.VINA_SEED, self.parameters.seed,
                     _EE.VINA_OUT, tmp_pdbqt_docked,
                     _EE.VINA_CENTER_X, str(search_space.center_x),
                     _EE.VINA_CENTER_Y, str(search_space.center_y),
                     _EE.VINA_CENTER_Z, str(search_space.center_z),
                     _EE.VINA_SIZE_X, str(search_space.size_x),
                     _EE.VINA_SIZE_Y, str(search_space.size_y),
                     _EE.VINA_SIZE_Z, str(search_space.size_z),
                     _EE.VINA_NUM_MODES, self.parameters.number_poses]

        execution_result = self._ADV_executor.execute(command=_EE.VINA,
                                                      arguments=arguments,
                                                      check=True)
        self._delay4file_system(path=tmp_pdbqt_docked)

        # translate the parsed output PDBQT into an SDF
        arguments = [tmp_pdbqt_docked,
                     _BEE.OBABLE_INPUTFORMAT_PDBQT,
                     _BEE.OBABEL_OUTPUT_FORMAT_SDF,
                     "".join([_BEE.OBABEL_O, output_path_sdf])]
        self._OpenBabel_executor.execute(command=_BEE.OBABEL,
                                         arguments=arguments,
                                         check=False)
        self._delay4file_system(path=output_path_sdf)

    def write_docked_ligands(self, path, mode="all"):
        """This method overrides the parent class, docker.py write_docked_ligands method. This method writes docked
        ligands binding poses and conformers to a file. There is the option to output the best predicted binding pose
        per ligand, the best predicted binding pose per enumeration, or all the predicted binding poses

        :param path: Contains information on results output path
        :type path: string
        :param mode: Determines whether the output contains the best predicted binding pose per ligand, the best
            predicted binding pose per enumeration, or all the predicted binding poses
        :type mode: string, optional, default value is "all". Other possible values are "best_per_ligand" and
            "best_per_enumeration"
        :raises DockingRunFailed Error: This error is raised if the docking run has not been performed
        :raises OpenEye (OE) Fatal Error: This error is raised if the output file was unable to be created. Issues may
            be due to problems with the ligand structure
        :raises ValueError: This error is raised if the ligands are neither RDkit nor OpenEye readable
        """
        self._write_docked_ligands(path, mode, mol_type=_LP.TYPE_RDKIT)
