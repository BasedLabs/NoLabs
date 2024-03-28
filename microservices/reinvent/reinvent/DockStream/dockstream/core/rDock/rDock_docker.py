import os
import tempfile
import shutil
import multiprocessing
from copy import deepcopy
from typing import Optional, List, Any

import rdkit.Chem as Chem
from pydantic import BaseModel
from typing_extensions import Literal

from dockstream.core.Schrodinger.Glide_docker import Parallelization
from dockstream.core.docker import Docker
from dockstream.core.rDock.rDock_result_parser import rDockResultParser
from dockstream.utils.enums.logging_enums import LoggingConfigEnum
from dockstream.utils.execute_external.rDock import rDockExecutor
from dockstream.utils.enums.rDock_enums import rDockExecutablesEnum, rDockDockingConfigurationEnum, rDockRbdockOutputEnum
from dockstream.utils.enums.RDkit_enums import RDkitLigandPreparationEnum
from dockstream.utils.general_utils import gen_temp_file

from dockstream.utils.translations.molecule_translator import MoleculeTranslator
from dockstream.utils.dockstream_exceptions import DockingRunFailed


_LE = LoggingConfigEnum()
_LP = RDkitLigandPreparationEnum()
_CE = rDockDockingConfigurationEnum()
_EE = rDockExecutablesEnum()
_ROE = rDockRbdockOutputEnum()


class rDockParameters(BaseModel):
    prefix_execution: Optional[str] = None
    binary_location: Optional[str] = None
    parallelization: Optional[Parallelization]
    rbdock_prm_paths: List[str]
    number_poses: int

    def get(self, key: str) -> Any:
        """Temporary method to support nested_get"""
        return self.dict()[key]

    
class rDock(Docker):
    """Interface to the "rDock" backend."""
    backend: Literal["rDock"] = "rDock"
    parameters: rDockParameters

    _rDock_executor: rDockExecutor = None

    class Config:
        underscore_attrs_are_private = True

    def __init__(self, **run_parameters):
        super().__init__(**run_parameters)

    def _initialize_executor(self):
        """Initialize the executor for all "rDock" related calls and also check if it is available."""
        if self._rDock_executor is None:
            self._rDock_executor = rDockExecutor(
                prefix_execution=self.parameters.prefix_execution,
                binary_location=self.parameters.binary_location
            )
        if not self._rDock_executor.is_available():
            raise DockingRunFailed("Cannot initialize rDock docker, as rDock backend is not available - abort.")
        self._rDock_executor.set_env_vars()
        self._logger.log(f"Checked rDock backend availability (prefix_execution={self.parameters.prefix_execution}).", _LE.DEBUG)

    def _get_score_from_conformer(self, conformer):
        return float(conformer.GetProp(_ROE.SCORE))

    def add_molecules(self, molecules: list):
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

            # generate temporary input file and output directory into which "rbdock" will deposit the poses
            cur_tmp_output_dir = tempfile.mkdtemp()
            cur_tmp_sdf = gen_temp_file(prefix=str(start_index), suffix=".sdf", dir=cur_tmp_output_dir)

            # write-out the temporary input file
            one_written = False
            writer = Chem.SDWriter(cur_tmp_sdf)
            for ligand in sublist:
                # initialize all ligands (as they could have failed)
                if ligand.get_molecule() is not None:
                    mol = deepcopy(ligand.get_molecule())
                    one_written = True
                    mol.SetProp("_Name", ligand.get_identifier())
                    writer.write(mol)
            writer.close()
            if one_written is False:
                if os.path.isdir(cur_tmp_output_dir):
                    shutil.rmtree(cur_tmp_output_dir)
                continue

            tmp_output_dirs.append(cur_tmp_output_dir)
            tmp_input_sdf_paths.append(cur_tmp_sdf)
            tmp_output_sdf_paths.append('.'.join([cur_tmp_output_dir, "sd"]))
        return tmp_output_dirs, tmp_input_sdf_paths, tmp_output_sdf_paths

    def _dock(self, number_cores):

        self._initialize_executor()

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

            # run in parallel
            processes = []
            for chunk_index in range(len(tmp_output_dirs)):
                p = multiprocessing.Process(target=self._dock_subjob, args=(tmp_input_sdf_paths[chunk_index],
                                                                            tmp_output_dirs[chunk_index],
                                                                            tmp_output_sdf_paths[chunk_index]))
                processes.append(p)
                p.start()
            for p in processes:
                p.join()

            # add the number of input sublists rather than the output temporary folders to account for cases where
            # entire sublists failed to produce an input structure
            sublists_submitted += len(cur_slice_sublists)

            # load the chunks and recombine the result; add conformations
            for chunk_index in range(len(tmp_output_dirs)):
                if not os.path.isfile(tmp_output_sdf_paths[chunk_index]) or os.path.getsize(tmp_output_sdf_paths[chunk_index]) == 0:
                    continue

                # do not sanitize, because rDock sometimes produces stuff that cannot be kekulized
                for molecule in Chem.SDMolSupplier(tmp_output_sdf_paths[chunk_index], sanitize=False, removeHs=False):
                    # it can happen, that ligands have "impossible chemistry" and will be loaded by RDkit as "None"
                    if molecule is None:
                        continue
                    cur_conformer_name = str(molecule.GetProp(_ROE.NAME))

                    # add molecule to the appropriate ligand
                    for ligand in self.ligands:
                        if ligand.get_identifier() == cur_conformer_name:
                            ligand.add_conformer(molecule)
                            break

            # clean-up
            for path in tmp_output_dirs:
                shutil.rmtree(path)
            self._log_docking_progress(number_done=sublists_submitted, number_total=number_sublists)

        # sort the conformers (best to worst), update their names to contain the conformer id and add tags
        # -> <ligand_number>:<enumeration>:<conformer_number>
        for ligand in self.ligands:
            ligand.set_conformers(sorted(ligand.get_conformers(),
                                         key=lambda x: float(x.GetProp(_ROE.SCORE)), reverse=False))
            ligand.add_tags_to_conformers()

        # log any docking fails
        self._docking_fail_check()

        # parse the result of the docking step
        result_parser = rDockResultParser([ligand.get_clone() for ligand in self.ligands])
        self._df_results = result_parser.as_dataframe()

        # docking flag
        self._docking_performed = True

    def _dock_subjob(self, input_path_sdf, output_dir_path, output_sdf_path):

        # set up arguments list and execute
        # for an explanation of the parameters, see "rDockExecutablesEnum"
        # TODO: support "ensemble docking" - currently, only the first entry is used
        arguments = [_EE.RBDOCK_R, self.parameters.rbdock_prm_paths[0],
                     _EE.RBDOCK_I, input_path_sdf,
                     _EE.RBDOCK_O, output_dir_path,
                     _EE.RBDOCK_N, str(self.parameters.number_poses),
                     _EE.RBDOCK_S, str(_EE.RBDOCK_S_DEFAULT),
                     _EE.RBDOCK_P, _EE.RBDOCK_P_DEFAULT]

        execution_result = self._rDock_executor.execute(command=_EE.RBDOCK,
                                                        arguments=arguments,
                                                        check=True)
        self._delay4file_system(path=output_sdf_path)
        self._logger.log(f"Finished sublist (input: {input_path_sdf}, output directory: {output_dir_path}).", _LE.DEBUG)

    def write_docked_ligands(self, path, mode="all"):
        self._write_docked_ligands(path, mode, mol_type=_LP.TYPE_RDKIT)
