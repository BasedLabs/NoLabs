import abc
import time
from copy import deepcopy
import multiprocessing
from enum import Enum
from typing import List, Optional, Union

import pandas as pd
from pydantic import BaseModel, PrivateAttr

from dockstream.loggers.docking_logger import DockingLogger
from dockstream.loggers.blank_logger import BlankLogger
from dockstream.utils.dockstream_exceptions import DockingRunFailed
from dockstream.utils.files_paths import generate_folder_structure

from dockstream.utils.parallelization.general_utils import split_into_sublists, get_progress_bar_string
from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum
from dockstream.utils.enums.docking_enum import DockingConfigurationEnum, ResultKeywordsEnum
from dockstream.utils.enums.logging_enums import LoggingConfigEnum
from dockstream.utils.general_utils import *

_DE = DockingConfigurationEnum()
_RK = ResultKeywordsEnum()
_LE = LoggingConfigEnum()
_LPE = LigandPreparationEnum()


class OutputMode(str, Enum):
    """Output mode.

    Determines whether the output contains the best predicted binding pose per ligand,
    the best predicted binding pose per enumeration,
    or all the predicted binding poses.
    """

    ALL = "all"
    BESTPERLIGAND = "best_per_ligand"
    BESTPERENUMERATION = "best_per_enumeration"


class Poses(BaseModel):
    poses_path: str
    overwrite: bool = False
    mode: OutputMode = OutputMode.BESTPERLIGAND

    class Config:
        use_enum_values = True


class Scores(BaseModel):
    scores_path: str
    overwrite: bool = False
    mode: OutputMode = OutputMode.BESTPERLIGAND

    class Config:
        use_enum_values = True


class Output(BaseModel):
    poses: Poses
    scores: Scores


class Docker(BaseModel, metaclass=abc.ABCMeta):
    """Virtual class implementing the interface to the actual docking backends."""

    input_pools: Union[str, List[str]]
    output: Optional[Output]
    run_id: Optional[str]

    ligands: List = []

    _logger = PrivateAttr()
    _logger_blank = PrivateAttr()

    _df_results = PrivateAttr()
    _run_parameters = PrivateAttr()
    _docking_performed = PrivateAttr()

    class Config:
        underscore_attrs_are_private = True
        extra = "allow"

    def __init__(self, **data):
        super().__init__(**data)

        self._logger = DockingLogger()
        self._logger_blank = BlankLogger()

        # the results pandas dataframe will be filled by the docking backends
        self._df_results = None

        # store the specific parameters for this very docking run for easy access later; the others are ignored
        self._run_parameters = data

        # prepare the list
        self.ligands = []

        self._docking_performed = False

    def add_molecules(self, molecules: list):
        """This method appends prepared ligands for docking to a list. It must be overrode by an add_molecules method
        in each backend (ex. Schrodinger Glide)
        :raises NotImplementedError: This error is raised if add_molecules is not implemented in a given backend
        """
        raise NotImplementedError

    def dock(self):
        """This method takes a given backend (ex. Schrodinger Glide) and docks the prepared ligands
        :raises Exception: An exception is raised if the ligands list is empty. The ligand list must first me
            populated using the backend's add_molecules method before docking can commence
        """
        if len(self.ligands) == 0:
            raise Exception("Add molecules to dock first.")

        # delete conformers
        for ligand in self.ligands:
            ligand.clear_conformers()

        # prepare the parallelization and set the number of cores to be used
        number_cores = nested_get(self._run_parameters, [_DE.PARAMS,
                                                         _DE.PARALLELIZATION,
                                                         _DE.PARALLELIZATION_NUMBER_CORES], default=1)
        if number_cores == 0:
            number_cores = 1
        elif number_cores < 0:
            # subtract the number of cores (neg. value, thus add up) from total number of cores, e.g. -1 will
            # use all available cores minus 1
            number_cores = multiprocessing.cpu_count() + number_cores

        # call the backend-specific, overloaded docking routine
        self._dock(number_cores=number_cores)

    def _dock(self, number_cores):
        raise NotImplementedError

    def has_result(self) -> bool:
        """This method returns whether the pandas dataframe has been populated with docking data by a given backend
        (ex. Schrodinger Glide).
        :return: boolean, either True or False
        """
        if self._df_results is not None:
            return True
        else:
            return False

    def get_result(self) -> pd.DataFrame:
        """This method returns the pandas dataframe results from a given docking run
        :return: pd.DataFrame, pandas dataframe containing docking results
        :raises DockingRunFailed Error: This error is raised if the docking run has yet to be performed
        """
        if self.has_result():
            return self._df_results
        else:
            raise DockingRunFailed("Execute method dock() before accessing the results.")

    def _log_docking_progress(self, number_done, number_total):
        self._logger.log(get_progress_bar_string(number_done, number_total, length=65),
                         _LE.INFO)

    def get_sublists_for_docking(self, number_cores, enforce_singletons=False):
        """This method splits ligands into sublists for docking to take advantage of parallel computing using >1
        processing core on your computer.
        :param number_cores: Number of cores you want to allocate to the docking job. Ensure at least 1 core
            is kept free so other tasks can still run on your computer (ex. if you have 8 cores, allocate at most 7)
        :type number_cores: int
        :param enforce_singletons: Determines whether a singleton pattern is to be enforced. This would split ligands
        into sublists such that each sublist only contains 1 ligand. This is necessary for AutoDock Vina
        :type enforce_singletons: boolean, default value is False. Possible values include True or False
        :return: split_into_sublists containing ligands split into sublists for subsequent parallel docking
        """
        # if every sublist should have exactly one member, split it (e.g. for AutoDock Vina)
        if enforce_singletons:
            return split_into_sublists(input_list=self.ligands, partitions=None, slice_size=1)

        # decide how to slice the ligand list depending on whether a maximum length is defined or not
        max_compounds_per_subjob = \
            nested_get(self._run_parameters, [_DE.PARAMS,
                                             _DE.PARALLELIZATION,
                                             _DE.PARALLELIZATION_MAXCOMPOUNDSPERSUBJOB],
                      default=0)
        if max_compounds_per_subjob > 0:
            slice_size = min(max_compounds_per_subjob,
                             len(self.ligands))
            return split_into_sublists(input_list=self.ligands, partitions=None, slice_size=slice_size)
        else:
            # split the ligands into as many cores as available
            partitions = min(number_cores, len(self.ligands))
            return split_into_sublists(input_list=self.ligands, partitions=partitions, slice_size=None)

    def get_docked_ligands(self):
        """This method returns a list of the docked ligand poses from a given docking run
        :raises DockingRunFailed Error: This error is raised if the docking has not been run yet
        :return: list containing all the docked ligand poses
        """
        return [ligand.get_clone() for ligand in self.ligands]

    def write_docked_ligands(self, path, mode="all"):
        """This method appends writes docked ligands binding poses and conformers to a file. There is the option
        to output the best predicted binding pose per ligand, the best predicted binding pose per enumeration,
        or all the predicted binding poses. This method must be overrode in some backends (ex. Schrodinger Glide)

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
        raise NotImplementedError

    def _is_best_per_enumeration(self, conformer, mol_type):
        # assume, that the conformers are already ordered and the last position represents the conformer ID
        if mol_type == _LPE.TYPE_RDKIT:
            if not conformer.GetProp("_Name").endswith(":0"):
                return False
        elif mol_type == _LPE.TYPE_OPENEYE:
            if not conformer.GetTitle().endswith(":0"):
                return False
        else:
            self._logger.log(f"Unknown molecule type {mol_type}.", _LE.ERROR)
            raise ValueError(f"Unknown molecule type {mol_type}.")
        return True

    def _get_ligand_id_from_conformer(self, mol_type: str, conformer):
        if mol_type == _LPE.TYPE_OPENEYE:
            return conformer.GetTitle().split(':')[0]
        else:
            return conformer.GetProp("_Name").split(':')[0]

    def _get_best_conformer_per_ligand(self, conformers: list, mol_type: str, ligand_ids: list) -> list:
        # group conformers by ligand IDs
        dict_grouped = {str(ligand_id): [] for ligand_id in ligand_ids}
        for conformer in conformers:
            cur_id = self._get_ligand_id_from_conformer(mol_type=mol_type, conformer=conformer)
            dict_grouped[str(cur_id)].append(conformer)

        # get the best conformer per ligand
        selected_conformers = []
        for ligand_id in ligand_ids:
            list_conf = dict_grouped[str(ligand_id)]
            if len(list_conf) > 0:
                list_conf = self._sort_conformers(conformers=list_conf)
                selected_conformers.append(list_conf[0])
        return selected_conformers

    def _sort_conformers(self, conformers: list, best="min") -> list:
        if best == "min":
            return sorted(conformers, key=lambda c: self._get_score_from_conformer(conformer=c))
        elif best == "max":
            return sorted(conformers, key=lambda c: self._get_score_from_conformer(conformer=c), reverse=True)
        else:
            self._logger.log(f"Parameter best must be either \"min\" or \"max\" (value {best} unknown).",
                             _LE.EXCEPTION)
            raise ValueError

    def _select_conformers(self, mode, mol_type):
        ligands = [deepcopy(lig) for lig in self.ligands]
        selected_conformers = []

        # extract all conformers for all ligands
        for ligand in ligands:
            for conformer in ligand.get_conformers():
                selected_conformers.append(conformer)
        if mode == _DE.OUTPUT_MODE_ALL:
            return selected_conformers

        # filter down to "best_per_enumeration"
        selected_conformers = [conf for conf in selected_conformers if self._is_best_per_enumeration(conformer=conf,
                                                                                                     mol_type=mol_type)]
        if mode == _DE.OUTPUT_MODE_BESTPERENUMERATION:
            return selected_conformers

        # filter down to "best_per_ligand"
        selected_conformers = self._get_best_conformer_per_ligand(conformers=selected_conformers,
                                                                  mol_type=mol_type,
                                                                  ligand_ids=list(set([ligand.get_ligand_number() for ligand in ligands])))
        return selected_conformers

    def _write_docked_ligands(self, path, mode, mol_type):
        if not self._docking_performed:
            raise DockingRunFailed("Do the docking first.")
        selected_conformers = self._select_conformers(mode=mode, mol_type=mol_type)

        # generate folder structure, if not available
        generate_folder_structure(filepath=path)

        if mol_type == _LPE.TYPE_RDKIT:
            import rdkit.Chem as Chem
            writer = Chem.SDWriter(path)
            for conformer in selected_conformers:
                writer.write(conformer)
            writer.close()
        elif mol_type == _LPE.TYPE_OPENEYE:
            import openeye.oechem as oechem
            ofs = oechem.oemolostream()
            ofs.SetFormat(oechem.OEFormat_SDF)
            if ofs.open(path):
                for conformer in selected_conformers:
                    oechem.OEWriteMolecule(ofs, conformer)
            else:
                oechem.OEThrow.Fatal("Unable to create specified output file.")
            ofs.close()
        self._logger.log(f"Wrote docked ligands to file {path}.", _LE.DEBUG)

    def _docking_fail_check(self):
        """This method checks for docking failures for two cases:
           1) if any ligand and all its enumerations fails to dock
           2) any enumeration fails to dock
           Relevant messages are logged"""
        self._logger.log(f"Attempted to dock {len(self.ligands)} molecules.", _LE.DEBUG)
        # check if there are cases where a ligand completely fails to dock (i.e. all its enumerations fail)
        ligand_num_tracker = []
        total_ligand_fails = 0

        for ligand in self.ligands:
            current_num = ligand.get_ligand_number()
            if current_num not in ligand_num_tracker:
                ligand_num_tracker.append(current_num)
                fail_bool = False
                for docked_ligand in self.ligands:
                    if docked_ligand.get_ligand_number() == current_num:
                        if len(docked_ligand.get_conformers()) == 0:
                            smiles = docked_ligand.get_original_smile()
                            fail_bool = True
                        else:
                            fail_bool = False
                            break

                if fail_bool:
                    # this block only runs if the ligand and all its enumerations failed to dock
                    total_ligand_fails += 1
                    self._logger.log(f"Ligand {ligand.get_ligand_number()} with SMILES: {smiles} and all its enumerations failed to dock.", _LE.DEBUG)

        self._logger.log(f"{total_ligand_fails} ligand(s) completely failed to dock", _LE.DEBUG)

        # keep track of the number of enumerated ligands which failed to dock
        not_docked = len([conf for conf in self.ligands if len(conf.get_conformers()) == 0])
        self._logger.log(f"{not_docked} ligand enumeration(s) failed to dock (did not return a pose and score)", _LE.DEBUG)

    def write_result(self, path, mode="all"):
        """This method writes the docking results to a csv file. There is the option to write out either the best
        predicted binding pose per enumeration or all the predicted binding poses. Output for the best predicted
        binding pose per ligand has yet to be implemented

        :param path: Contains information on results output path
        :type path: string
        :param mode: Determines whether the output contains the best predicted binding pose per ligand, the best
            predicted binding pose per enumeration, or all the predicted binding poses
        :type mode: string, optional, default value is "all". Other possible value is "best_per_enumeration"
        """
        return self._write_result(path=path, mode=mode)

    def _write_result(self, path, mode="all", best="min"):
        if self._df_results is not None:
            if self._df_results.empty:
                self._logger.log("Generated dataframe is empty, skipping write-out (this probably means all poses were rejected).",
                                 _LE.WARNING)
            else:
                if mode == _DE.OUTPUT_MODE_ALL:
                    df_buffer = self._df_results.copy()
                elif mode == _DE.OUTPUT_MODE_BESTPERENUMERATION:
                    df_buffer = self._df_results.copy()
                    df_buffer = df_buffer[df_buffer[_RK.DF_LOWEST_CONFORMER]]
                elif mode == _DE.OUTPUT_MODE_BESTPERLIGAND:
                    df_buffer = self._df_results.copy()
                    df_buffer = df_buffer[df_buffer[_RK.DF_LOWEST_CONFORMER]]
                    if best == "min":
                        df_buffer = df_buffer.loc[df_buffer.groupby(_RK.DF_LIGAND_NUMBER)[_RK.DF_SCORE].idxmin()]
                    elif best == "max":
                        df_buffer = df_buffer.loc[df_buffer.groupby(_RK.DF_LIGAND_NUMBER)[_RK.DF_SCORE].idxmax()]
                    else:
                        self._logger.log(f"Parameter best must be either \"min\" or \"max\" (value {best} unknown).",
                                         _LE.EXCEPTION)
                        raise ValueError
                else:
                    self._logger.log(f"Score output mode \"{mode}\" is unknown - write-out of scores failed.",
                                     _LE.ERROR)
                    raise DockingRunFailed()

                # generate folder structure, if not available
                generate_folder_structure(filepath=path)

                df_buffer.to_csv(path_or_buf=path,
                                 sep=',',
                                 na_rep='',
                                 header=True,
                                 index=False,
                                 mode='w',
                                 quoting=None)
                self._logger.log(f"Wrote result of docking run to file {path} with mode set to \"{mode}\" ({df_buffer.shape[0]} rows).",
                                 _LE.DEBUG)

    def get_scores(self, best_only):
        """This method returns a list containing all docking scores. This method allows returning the best docking
        scores only. "Best" can mean the minimum or maximum values for this given scoring function. By default,
        it will return the minimum values (greatest predicted binding energy). Returning all the docking scores
        (of different poses for instance) is also possible if best only is not enforced. This method must be overrode
        in some backends (ex. GOLD)

        :param best_only: Determines whether the best (either minimum or maximum) docking scores are returned
        :type best_only: boolean, True of False
        :return: list of returned docking scores
        :raises ValueError: If best_only is True but neither "min" nor "max" was specified, a ValueError is raised
        """
        return self._get_scores(best_only, "min")

    def _get_scores(self, best_only, best):
        # note, that here one needs to combine all Ligand objects with the same ligand_number together (they should
        # differ in their enumeration value); also note that if multiple conformers are available, they have been
        # sorted already in descending performance
        if not self._docking_performed:
            raise DockingRunFailed("Do the docking first.")

        # combine scores of all enumerations (Ligand objects) into lists
        ligand_numbers = list(set([ligand.get_ligand_number() for ligand in self.ligands]))
        buffer_list = []
        for ligand_number in ligand_numbers:
            cur_ligand_list = []
            for ligand in self.ligands:
                if ligand_number != ligand.get_ligand_number():
                    continue
                for conformer in ligand.get_conformers():
                    cur_ligand_list.append(self._get_score_from_conformer(conformer))
            buffer_list.append(cur_ligand_list)

        # empty list -> no valid docking, return "NA"
        # otherwise store either the best score per list or all values
        result_list = []
        for ligand_scores_list in buffer_list:
            if len(ligand_scores_list) == 0:
                result_list.append(_RK.FIXED_VALUE_NA)
            else:
                if best_only:
                    if best == "min":
                        result_list.append(min(ligand_scores_list))
                    elif best == "max":
                        result_list.append(max(ligand_scores_list))
                    else:
                        self._logger.log(f"Parameter best must be either \"min\" or \"max\" (value {best} unknown).",
                                         _LE.EXCEPTION)
                        raise ValueError
                else:
                    result_list = result_list + ligand_scores_list
        return result_list

    def _get_score_from_conformer(self, conformer):
        raise NotImplementedError

    @staticmethod
    def apply_prefix_to_filename(path: str, output_prefix: str):
        """This method modifies the output path directory by concatenating an output_prefix at the front

        :param path:: The base output path pre-modification
        :type path: string
        :param output_prefix: prefix to be concatenated at the front to modify the output path
        :type output_prefix: string
        :return: string path with output_prefix concatenated at the front if path is not None.
            Otherwise, return the base path
        """
        if output_prefix is not None:
            return os.path.join(os.path.dirname(path), "".join([output_prefix, os.path.basename(path)]))
        else:
            return path

    @staticmethod
    def update_path_to_unused(path: str) -> str:
        """This function takes a path, checks if it is already in use and if so generates a higher-indexed version
        of the path (e.g. "<the_path>/output.sdf" will become "<the_path>/0001_output.sdf" and so on).

        :param path: The base output path
        :type path: string
        :return: string, output path
        """

        # check if file specified by path exists and if not keep the input path
        # note: there seems to be a Python bug that sometimes returns False for "path.isfile" even though the file
        #       exists; use "exists" instead and keep an eye on that
        if not os.path.exists(path):
            return path

        # file exists; find next path that is suitable
        index = 0
        filename = os.path.basename(path)
        while index <= 9999:
            index = index + 1
            prefix = f"{index:04}"
            prosp_path = os.path.join(os.path.dirname(path), '_'.join([prefix, filename]))
            if os.path.exists(prosp_path):
                continue
            else:
                return prosp_path
        raise DockingRunFailed("Could not find path replacement.")

    def _wait_until_file_generation(self, path, interval_sec=1, maximum_sec=None) -> bool:
        counter = 0
        while not os.path.exists(path):
            # wait for an interval
            time.sleep(interval_sec)
            counter = counter + 1

            # if there's time left, proceed
            if maximum_sec is not None and counter * interval_sec >= maximum_sec:
                break
        if os.path.exists(path):
            return True
        else:
            return False

    def _delay4file_system(self, path) -> bool:
        return self._wait_until_file_generation(path=path, interval_sec=1, maximum_sec=10)
