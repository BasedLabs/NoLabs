import os
import shutil
import multiprocessing
import tempfile
from rdkit import Chem

from typing import Optional
from typing_extensions import Literal   # Required for Python 3.7. From 3.8 Literal is in typing.
from pydantic import BaseModel, PrivateAttr

from dockstream.core.RDkit.RDkit_ligand_preparator import RDkitLigandPreparator

from copy import deepcopy
from dockstream.core.ligand.ligand import Ligand, get_enumerations_for_ligand

from dockstream.core.ligand_preparator import LigandPreparator, _LE
from dockstream.utils.general_utils import gen_temp_file
from dockstream.utils.parallelization.general_utils import split_into_sublists, get_progress_bar_string

from dockstream.loggers.blank_logger import BlankLogger

from dockstream.utils.dockstream_exceptions import LigandPreparationFailed

from dockstream.utils.translations.molecule_translator import MoleculeTranslator
from dockstream.utils.smiles import to_mol, to_smiles

from dockstream.utils.execute_external.Omega import OmegaExecutor
from dockstream.utils.enums.Omega_enums import OmegaExecutablesEnum, OmegaOutputEnum

_OE = OmegaExecutablesEnum()
_OO = OmegaOutputEnum()


class Parallelization(BaseModel):
    number_cores: int = 1
    max_compounds_per_subjob: Optional[int] = None


class OmegaLigandPreparatorParameters(BaseModel):
    prefix_execution: Optional[str] = None
    binary_location: Optional[str] = None
    parallelization: Optional[Parallelization] = Parallelization()
    mode: str = "classic"


class OmegaLigandPreparator(LigandPreparator, BaseModel):
    """Class that acts as an interface to the "OMEGA" executable to prepare ligands."""

    type: Literal["Omega"] = "Omega"
    parameters: OmegaLigandPreparatorParameters = OmegaLigandPreparatorParameters()

    class Config:
        underscore_attrs_are_private = True

    _omega_executor: OmegaExecutor = PrivateAttr()
    _logger_blank = PrivateAttr()

    def __init__(self, **data):
        super().__init__(**data)

        self._logger_blank = BlankLogger()

        self._omega_executor = OmegaExecutor(prefix_execution=self.parameters.prefix_execution,
                                             binary_location=self.parameters.binary_location)
        if not self._omega_executor.is_available():
            raise LigandPreparationFailed("Cannot initialize OMEGA backend - abort.")
        self._logger.log(f"Checked OMEGA backend availability (prefix_execution={self.parameters.prefix_execution}, binary_location={self.parameters.binary_location}).",
                         _LE.DEBUG)

    def _initialize_ligands(self):
        super()._initialize_ligands()

        # here, it is guaranteed that the molecules are given as list of "Ligand" objects
        mol_type = self.ligands[0].get_mol_type()
        if mol_type is None:
            # no molecule attached -> transform smiles into rdkit molecules
            for lig in self.ligands:
                mol = to_mol(lig.get_smile())
                lig.set_molecule(mol)
                lig.set_mol_type(_OE.TYPE_OMEGA)
        else:
            mol_trans = MoleculeTranslator(self.ligands)
            self.ligands = mol_trans.get_as_rdkit()
        self._logger.log(f"Stored {len(self.ligands)} OMEGA molecules in ligands (RDkit type).",
                         _LE.DEBUG)

    def _load_references(self):
        references = []
        ref_format = self.align.reference_format.upper()
        for path in self.align.reference_paths:
            if ref_format == _OE.ALIGN_REFERENCE_FORMAT_PDB:
                ref_mol = Chem.MolFromPDBFile(path, sanitize=True)
                ref_mol.SetProp("_Name", os.path.basename(path))
                references.append(ref_mol)
            elif ref_format == _OE.ALIGN_REFERENCE_FORMAT_SDF:
                mol_supplier = Chem.SDMolSupplier(path)
                for mol in mol_supplier:
                    references.append(mol)
            else:
                raise IOError("Specified format not supported.")
        if len(references) == 0:
            raise LigandPreparationFailed("No reference molecules could be loaded with path(s) specified.")
        self._references = references
        self._logger.log(f"Stored {len(references)} reference molecules.", _LE.DEBUG)

    def _get_RDkit_aligner(self, conf, ligands):
        return RDkitLigandPreparator(ligands=ligands, **conf)

    def _prepare_omega_arguments(self) -> list:
        arguments_list = []

        # add the OMEGA mode specified
        if self.parameters.mode not in [_OE.CLASSIC, _OE.MACROCYCLE, _OE.ROCS, _OE.POSE, _OE.DENSE]:
            raise ValueError(f"{self.parameters.mode} mode is not valid. Supported modes are:"
                             f"'classic', 'macrocycle, 'rocs', 'pose', 'dense'")
        else:
            arguments_list.append(self.parameters.mode)

        # Add other parameters to expose below

        return arguments_list

    def _get_number_cores(self):
        # prepare the parallelization and set the number of cores to be used
        number_cores = self.parameters.parallelization.number_cores
        if number_cores == 0:
            number_cores = 1
        elif number_cores < 0:
            # subtract the number of cores (neg. value, thus add up) from total number of cores, e.g. -1 will
            # use all available cores minus 1
            number_cores = multiprocessing.cpu_count() + number_cores
        return number_cores

    def _get_sublists_for_embedding(self, number_cores, enforce_singletons=False):
        if enforce_singletons:
            return split_into_sublists(input_list=self.ligands, partitions=None, slice_size=1)

        # decide how to slice the ligand list depending on whether a maximum length is defined or not
        if self.parameters.parallelization.max_compounds_per_subjob is not None:
            slice_size = min(max(self.parameters.parallelization.max_compounds_per_subjob, 1),
                             len(self.ligands))
            return split_into_sublists(input_list=self.ligands, partitions=None, slice_size=slice_size)
        else:
            # split the ligands into as many cores as available
            partitions = min(number_cores, len(self.ligands))
            return split_into_sublists(input_list=self.ligands, partitions=partitions, slice_size=None)

    def _generate_temporary_input_output_files(self, start_indices, sublists):
        # in case singletons are handed over, wrap them in a list for "zipping" later
        if not isinstance(start_indices, list):
            start_indices = [start_indices]
        if not isinstance(sublists, list):
            sublists = [sublists]

        tmp_output_dirs = []
        tmp_input_smi_paths = []
        tmp_output_sdf_paths = []
        for start_index, sublist in zip(start_indices, sublists):
            # generate temporary input files and output directory
            cur_tmp_output_dir = tempfile.mkdtemp()
            tmp_output_dirs.append(cur_tmp_output_dir)
            cur_tmp_sdf = gen_temp_file(prefix=str(start_index), suffix=".sdf", dir=cur_tmp_output_dir)
            cur_tmp_smi = gen_temp_file(prefix=str(start_index), suffix=".smi", dir=cur_tmp_output_dir)
            tmp_input_smi_paths.append(cur_tmp_smi)

            # write smiles to temporary file as "OMEGA" backend
            with open(cur_tmp_smi, 'w') as f:
                for lig in sublist:
                    f.write(lig.get_smile() + " " + lig.get_identifier() + "\n")

            # add the path to which "_dock_subjob()" will write the result SDF
            output_sdf_path = gen_temp_file(prefix=str(start_index), suffix="_result.sdf", dir=cur_tmp_output_dir)
            tmp_output_sdf_paths.append(output_sdf_path)

        return tmp_output_dirs, tmp_input_smi_paths, tmp_output_sdf_paths

    def _log_docking_progress(self, number_done, number_total):
        self._logger.log(get_progress_bar_string(number_done, number_total, length=65), _LE.INFO)

    def generate3Dcoordinates(self):
        """Method to generate 3D coordinates, in case the molecules have to be built from SMILES."""

        number_cores = self._get_number_cores()
        start_indices, sublists = self._get_sublists_for_embedding(number_cores=number_cores)
        number_sublists = len(sublists)
        self._logger.log(f"Split ligands into {number_sublists} sublists for embedding.", _LE.DEBUG)

        sublists_submitted = 0
        slices_per_iteration = min(number_cores, number_sublists)
        if isinstance(self.ligands[0].get_molecule(), Chem.Mol):
            while sublists_submitted < len(sublists):
                upper_bound_slice = min((sublists_submitted + slices_per_iteration), len(sublists))
                cur_slice_start_indices = start_indices[sublists_submitted:upper_bound_slice]
                cur_slice_sublists = sublists[sublists_submitted:upper_bound_slice]

                # generate paths and initialize molecules (so that if they fail, this can be covered)
                tmp_output_dirs, tmp_input_smi_paths, \
                tmp_output_sdf_paths = self._generate_temporary_input_output_files(cur_slice_start_indices, cur_slice_sublists)

                # run in parallel; wait for all subjobs to finish before proceeding
                processes = []
                for chunk_index in range(len(tmp_output_dirs)):
                    p = multiprocessing.Process(target=self._run_embedding_subjob,
                                                args=(tmp_input_smi_paths[chunk_index],
                                                      tmp_output_sdf_paths[chunk_index],
                                                      tmp_output_dirs[chunk_index]))
                    processes.append(p)
                    p.start()
                for p in processes:
                    p.join()

                # add the number of input sublists rather than the output temporary folders to account for cases where
                # entire sublists failed to produce an input structure
                sublists_submitted += len(cur_slice_sublists)

                # load and store the conformers; name it sequentially
                # note, that some backends require the H-coordinates (such as Glide) - so keep them!
                ligands_embedded = []
                first_subjob, first_mol, idx = True, False, 0
                for path_sdf_results in tmp_output_sdf_paths:
                    if not os.path.isfile(path_sdf_results):
                        continue
                    if os.path.getsize(path_sdf_results) == 0:
                        self._logger.log(f"skipped output file as it is empty - typically, this indicates errors during OMEGA embedding.",
                                         _LE.DEBUG)
                        continue

                    mol_supplier = Chem.SDMolSupplier(path_sdf_results, removeHs=False)

                    # below boolean trackers are used to ensure only the first conformer generated per ligand is carried forward
                    # when parallelizing. The existence of many subjobs necessitate the safeguards below to ensure proper handling
                    if first_subjob:
                        first_subjob = False
                        idx = 0
                    else:
                        first_mol = True

                    for mol in mol_supplier:
                        if mol is not None and mol.HasProp("_Name"):

                            lig_id, enum_id = mol.GetProp("_Name").split(":")

                            if idx == 0 or first_mol:
                                first_mol = False
                                ligands_embedded.append(
                                    Ligand(smile="", ligand_number=int(lig_id), enumeration=int(enum_id),
                                           molecule=mol, mol_type=_OE.TYPE_OMEGA))
                                idx += 1

                            elif ligands_embedded[idx-1].get_ligand_number() == int(lig_id):
                                continue

                            else:
                                ligands_embedded.append(
                                    Ligand(smile="", ligand_number=int(lig_id), enumeration=int(enum_id),
                                           molecule=mol, mol_type=_OE.TYPE_OMEGA))
                                idx += 1

                        else:
                            self._logger.log(f"Skipped molecule when loading as _Name property could not be found - typically, this indicates that OMEGA could not embed the molecule.", _LE.WARNING)
                # update internal (self.ligands) list of ligands with new molecules
                # TODO: check expand enumerations functionality for OMEGA
                self._expand_enumerations(ligands_embedded)

                # remove temporary files
                for path in tmp_output_dirs:
                    shutil.rmtree(path)
                self._log_docking_progress(number_done=sublists_submitted, number_total=number_sublists)

        # check success and failure with embedding
        failed = 0
        for ligand in self.ligands:
            if ligand.get_molecule() is None:
                failed += 1
                self._logger.log(f"Enumeration {ligand.get_identifier()} could not be embedded (smile: {ligand.get_smile()}).", _LE.DEBUG)
        if failed > 0:
            self._logger.log(f"It appears OMEGA could not embed all {len(self.ligands)} ligands ({failed} were not found) - compounds not embedded might be ignored at the docking stage.",_LE.WARNING)

    def _run_embedding_subjob(self, smi_ligand_path, path_sdf_results, tmp_output_dir):
        # 1) prepare "OMEGA" arguments
        arguments = self._prepare_omega_arguments()
        arguments = arguments + [_OE.IN, smi_ligand_path]
        arguments = arguments + [_OE.OUT, path_sdf_results]

        # 2) run "OMEGA" backend and add log file to "debug" mode logging
        result = self._omega_executor.execute(command=_OE.OMEGA,
                                              arguments=arguments,
                                              location=tmp_output_dir,
                                              check=False)
        self._logger.log(f"Executed OMEGA backend (output file: {path_sdf_results}).", _LE.DEBUG)
        for file in os.listdir(tmp_output_dir):
            if file.endswith(".fail"):
                with open(os.path.join(tmp_output_dir, file), "r") as f:
                    for line in f.readlines():
                        error_components = line.split(" ")
                        self._logger.log(f"It appears OMEGA failed to embed: {error_components[0]}, with ligand ID: {error_components[1]} "
                                         f", due to the following error: {' '.join(error_components[2:])}", _LE.DEBUG)

        # below is code to print the OMEGA log file into the internal outputted log file. It is commented
        # out at the moment due to its sheer size
        #for file in os.listdir(tmp_output_dir):
            #if file.endswith(".log"):
                #path_tmp_log = os.path.join(tmp_output_dir, file)
                #self._print_log_file(path=path_tmp_log)

    def _expand_enumerations(self, ligands_embedded):
        # store the generated conformations with the original ligands; if embedding failed, keep the old (emtpy) one
        new_ligand_list = []
        for ligand in self.ligands:
            enums = get_enumerations_for_ligand(ligands=ligands_embedded, ligand_id=ligand.get_ligand_number())
            if len(enums) == 0:
                new_ligand_list.append(ligand)
            else:
                for enum_id, enum in enumerate(enums):
                    buf_lig = deepcopy(ligand)
                    buf_lig.set_enumeration(enum_id)
                    buf_lig.set_molecule(enum.get_molecule())
                    buf_lig.set_mol_type(enum.get_mol_type())
                    buf_lig.set_smile(to_smiles(enum.get_molecule()))
                    new_ligand_list.append(buf_lig)
        self.ligands = new_ligand_list

    def _print_log_file(self, path):
        if os.path.isfile(path):
            with open(path, 'r') as log_file:
                log_file_raw = log_file.readlines()
                self._logger.log(f"Printing log file {path}:\n", _LE.DEBUG)
                for line in log_file_raw:
                    self._logger_blank.log(line.rstrip("\n"), _LE.DEBUG)
                self._logger_blank.log("", _LE.DEBUG)
                self._logger.log("--- End file", _LE.DEBUG)

    def align_ligands(self):
        self.ligands = self._align_ligands_with_RDkit_preparator(self.ligands)
