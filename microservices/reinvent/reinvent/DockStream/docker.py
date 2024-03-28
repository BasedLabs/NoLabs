#!/usr/bin/env python
#  coding=utf-8

import os
import sys
import warnings
import argparse

from dockstream.containers.docking_container import DockingContainer

from dockstream.core.rDock.rDock_docker import rDock
from dockstream.core.OpenEye.OpenEye_docker import OpenEye
from dockstream.core.Schrodinger.Glide_docker import Glide
from dockstream.core.AutodockVina.AutodockVina_docker import AutodockVina
from dockstream.core.OpenEyeHybrid.OpenEyeHybrid_docker import OpenEyeHybrid

from dockstream.utils.entry_point_functions.header import initialize_logging, set_environment
from dockstream.utils.entry_point_functions.embedding import embed_ligands
from dockstream.utils.entry_point_functions.write_out import handle_poses_writeout, handle_score_printing, \
                                                         handle_scores_writeout

from dockstream.utils.enums.docking_enum import DockingConfigurationEnum
from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum
from dockstream.utils.enums.logging_enums import LoggingConfigEnum

from dockstream.utils.files_paths import attach_root_path
from dockstream.utils.argparse_bool_extension import str2bool
from dockstream.utils.dockstream_exceptions import *


if __name__ == "__main__":

    # enums
    _LE = LoggingConfigEnum()
    _LP = LigandPreparationEnum()
    _DE = DockingConfigurationEnum()

    # get the input parameters and parse them
    parser = argparse.ArgumentParser(description="Implements entry point for the docking using one or multiple backends.")
    parser.add_argument("-conf", type=str, default=None, help="A path to an docking configuration file (JSON dictionary) that is to be executed.", required=True)
    parser.add_argument("-validation", type=str2bool, default=True, help="If set to False, this flag will prohibit a JSON Schema validation.")
    parser.add_argument("-silent", type=str2bool, default=False, help="If set, the program will silently execute without printing status updates.")
    parser.add_argument("-smiles", default=None, help="Use this flag to hand over the input SMILES over the command-line, separated by ';'.", type=str, required=False)
    parser.add_argument("-print_scores", action="store_true", help="Set this flag to activate linewise print-outs of the scores to the shell.")
    parser.add_argument("-print_all", action="store_true", help="Set this flag (together with \"-print_scores\") to print out the scores for all conformers, not just the best one.")
    parser.add_argument("-debug", action="store_true", help="Set this flag to activate the inbuilt debug logging mode (this will overwrite parameter \"-log_conf\", if set).")
    parser.add_argument("-log_conf", type=str, default=None, help="Set absolute path to a logger configuration other than the default stored at \"config/logging/default.json\".")
    parser.add_argument("-output_prefix", type=str, default=None, help="If specified, this prefix will be added to all output file names.")
    parser.add_argument("-input_csv", type=str, default=None, help="If set (a path to a CSV file), this will overwrite any input file specification in the configuration.")
    parser.add_argument("-input_csv_smiles_column", type=str, default=None, help="If \"-input_csv\" is set, you need to specify the column name with the smiles as well.")
    parser.add_argument("-input_csv_names_column", type=str, default=None, help="Optional name of the name column, if \"-input_csv\" is specified.")
    args, args_unk = parser.parse_known_args()

    if args.conf is None or not os.path.isfile(args.conf):
        raise Exception("Parameter \"-conf\" must be a relative or absolute path to a configuration JSON file.")
    if args.print_scores is False and args.print_all:
        raise Exception("Flag \"-print_scores\" must be activated in order to use \"-print_all\", see help message.")

    # set the logging configuration according to parameters
    if args.log_conf is None:
        args.log_conf = attach_root_path(_LE.PATH_CONFIG_DEFAULT)
    if args.debug:
        args.log_conf = attach_root_path(_LE.PATH_CONFIG_DEBUG)

    # initialize the docking Enum and get the configuration
    try:
        config = DockingContainer(conf=args.conf, validation=args.validation)
    except Exception as e:
        raise DockingRunFailed() from e

    # header: process the header once before actually executing anything
    # ---------
    logger = initialize_logging(config=config, task=_DE.DOCKING, _task_enum=_DE, log_conf_path=args.log_conf)
    set_environment(config=config, task=_DE.DOCKING, _task_enum=_DE, logger=logger)

    # check, if there are unknown arguments
    if len(args_unk) > 0:
        logger.log(f"Unknown arguments: {args_unk}.", _LE.WARNING)

    # anything related to Gold (CCDC) fails if the proper environment has not been loaded; now, load the modules here
    try:
        with warnings.catch_warnings(record=True) as w:
            from dockstream.core.Gold.Gold_docker import Gold
            if len(w) > 0:
                logger.log("Could not load CCDC / Gold docker - if another backend is being used, you can safely ignore this warning.", _LE.DEBUG)
    except Exception as e:
        logger.log(f"Could not load CCDC / Gold docker - if another backend is being used, you can safely ignore this warning. The exception message reads: {get_exception_message(e)}", _LE.WARNING)

    # ligand preparation: transform SMILES into embedded (and potentially aligned) molecules
    #                     note, that this step is in principle independent from the actual docking
    # ---------
    dict_pools = {}
    if _LP.LIGAND_PREPARATION in config[_DE.DOCKING].keys():

        # If single element (from GUI), wrap in a list.
        if not isinstance(config[_DE.DOCKING][_LP.LIGAND_PREPARATION][_LP.EMBEDDING_POOLS], list):
            config[_DE.DOCKING][_LP.LIGAND_PREPARATION][_LP.EMBEDDING_POOLS] = [
                config[_DE.DOCKING][_LP.LIGAND_PREPARATION][_LP.EMBEDDING_POOLS]
            ]

        # check, if input specification for the pools is to be overwritten from the command-line
        if args.input_csv is not None:
            if args.input_csv_smiles_column is None:
                raise ValueError("When using \"-input_csv\", you need to also specify \"-input_csv_smiles_column\".")

            pools_list = config[_DE.DOCKING][_LP.LIGAND_PREPARATION][_LP.EMBEDDING_POOLS]
            new_input = {_LP.INPUT_TYPE: _LP.INPUT_TYPE_CSV,
                         _LP.INPUT_PATH: args.input_csv,
                         _LP.INPUT_CSV_DELIMITER: _LP.INPUT_CSV_DELIMITER_DEFAULT,
                         _LP.INPUT_CSV_COLUMNS: {
                             _LP.INPUT_CSV_COLNAME_SMILES: args.input_csv_smiles_column
                         }}
            if args.input_csv_names_column is not None:
                new_input[_LP.INPUT_CSV_COLUMNS][_LP.INPUT_CSV_COLNAME_NAMES] = args.input_csv_names_column

            for pool in pools_list:
                pool[_LP.INPUT] = new_input

        # ligand preparation is to be performed
        for pool_number, pool in enumerate(config[_DE.DOCKING][_LP.LIGAND_PREPARATION][_LP.EMBEDDING_POOLS]):
            logger.log(f"Starting generation of pool {pool[_LP.POOLID]}.", _LE.INFO)
            try:
                prep = embed_ligands(smiles=args.smiles,
                                     pool_number=pool_number,
                                     pool=pool,
                                     logger=logger,
                                     ligand_number_start=0)
                dict_pools[pool[_LP.POOLID]] = prep.get_ligands()
            except Exception as e:
                logger.log(f"Failed in constructing pool {pool[_LP.POOLID]}.", _LE.EXCEPTION)
                logger.log(f"Exception reads: {get_exception_message(e)}.", _LE.EXCEPTION)
                raise LigandPreparationFailed
            else:
                logger.log(f"Completed construction of pool {pool[_LP.POOLID]}.", _LE.INFO)

    # docking: this is the actual docking step; ligands can be provided by the preparation step specified before or
    #          loaded from files
    # ---------
    dict_docking_runs = {}
    if _DE.DOCKING_RUNS in config[_DE.DOCKING].keys():

        # If single element (from GUI), wrap in a list.
        if not isinstance(config[_DE.DOCKING][_DE.DOCKING_RUNS], list):
            config[_DE.DOCKING][_DE.DOCKING_RUNS] = [config[_DE.DOCKING][_DE.DOCKING_RUNS]]

        # execute the docking runs specified
        for docking_run_number, docking_run in enumerate(config[_DE.DOCKING][_DE.DOCKING_RUNS]):
            logger.log(f"Starting docking run {docking_run[_DE.RUN_ID]}.", _LE.INFO)
            try:
                if docking_run[_DE.BACKEND] == _DE.BACKEND_RDOCK:
                    docker = rDock(**docking_run)
                elif docking_run[_DE.BACKEND] == _DE.BACKEND_OPENEYE:
                    docker = OpenEye(**docking_run)
                elif docking_run[_DE.BACKEND] == _DE.BACKEND_OPENEYEHYBRID:
                    docker = OpenEyeHybrid(**docking_run)
                elif docking_run[_DE.BACKEND] == _DE.BACKEND_GLIDE:
                    docker = Glide(**docking_run)
                elif docking_run[_DE.BACKEND] == _DE.BACKEND_GOLD:
                    docker = Gold(**docking_run)
                elif docking_run[_DE.BACKEND] == _DE.BACKEND_AUTODOCKVINA:
                    docker = AutodockVina(**docking_run)
                else:
                    raise Exception("Backend is unknown.")

                # merge all specified pools for this run together
                if isinstance(docking_run[_DE.INPUT_POOLS], str):
                    docking_run[_DE.INPUT_POOLS] = [docking_run[_DE.INPUT_POOLS]]
                for pool_id in docking_run[_DE.INPUT_POOLS]:
                    cur_pool = [lig.get_clone() for lig in dict_pools.get(pool_id)]
                    if cur_pool is None or len(cur_pool) == 0:
                        raise Exception("Could not find pool id during docking run or pool was empty.")
                    docker.add_molecules(molecules=cur_pool)

                # do the docking
                docker.dock()

                # if specified, save the poses and the scores and print the scores to "stdout"
                handle_poses_writeout(docking_run=docking_run, docker=docker, output_prefix=args.output_prefix)
                handle_scores_writeout(docking_run=docking_run, docker=docker, output_prefix=args.output_prefix)
                handle_score_printing(print_scores=args.print_scores,
                                      print_all=args.print_all,
                                      docker=docker,
                                      logger=logger)
            except Exception as e:
                logger.log(f"Failed when executing run {docking_run[_DE.RUN_ID]}.", _LE.EXCEPTION)
                logger.log(f"Exception reads: {get_exception_message(e)}.", _LE.EXCEPTION)
                raise DockingRunFailed() from e
            else:
                logger.log(f"Completed docking run {docking_run[_DE.RUN_ID]}.", _LE.INFO)

    sys.exit(0)
