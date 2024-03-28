#!/usr/bin/env python
#  coding=utf-8

import os
import sys
import warnings
import argparse
import tempfile
from shutil import copyfile

from dockstream.core.pdb_preparator import PDBPreparator
from dockstream.utils.dockstream_exceptions import *

from dockstream.core.rDock.rDock_target_preparator import rDockTargetPreparator
from dockstream.core.OpenEye.OpenEye_target_preparator import OpenEyeTargetPreparator
from dockstream.core.AutodockVina.AutodockVina_target_preparator import AutodockVinaTargetPreparator
from dockstream.containers.target_preparation_container import TargetPreparationContainer

from dockstream.utils.entry_point_functions.header import initialize_logging, set_environment

from dockstream.utils.enums.target_preparation_enum import TargetPreparationEnum
from dockstream.utils.enums.OpenEye_enums import OpenEyeTargetPreparationEnum
from dockstream.utils.enums.Gold_enums import GoldTargetPreparationEnum
from dockstream.utils.enums.AutodockVina_enums import AutodockTargetPreparationEnum
from dockstream.utils.enums.logging_enums import LoggingConfigEnum

from dockstream.utils.files_paths import attach_root_path
from dockstream.utils.argparse_bool_extension import str2bool

from dockstream.utils.general_utils import *


if __name__ == "__main__":

    # enums
    _LE = LoggingConfigEnum()
    _TP = TargetPreparationEnum()

    # get the input parameters and parse them
    parser = argparse.ArgumentParser(description="Implements entry point for the target preparation for one or multiple backends.")
    parser.add_argument("-conf", type=str, default=None, help="A path to an preparation configuration file (JSON dictionary) that is to be executed.")
    parser.add_argument("-validation", type=str2bool, default=True, help="If set to False, this flag will prohibit a JSON Schema validation.")
    parser.add_argument("-silent", type=str2bool, default=False, help="If set, the program will silently execute without printing status updates.")
    parser.add_argument("-debug", action="store_true", help="Set this flag to activate the inbuilt debug logging mode (this will overwrite parameter \"-log_conf\", if set).")
    parser.add_argument("-log_conf", type=str, default=None, help="Set absolute path to a logger configuration other than the default stored at \"config/logging/default.json\".")
    args = parser.parse_args()

    if args.conf is None or not os.path.isfile(args.conf):
        raise Exception("Parameter \"-conf\" must be a relative or absolute path to a configuration JSON file.")

    # set the logging configuration according to parameters
    if args.log_conf is None:
        args.log_conf = attach_root_path(_LE.PATH_CONFIG_DEFAULT)
    if args.debug:
        args.log_conf = attach_root_path(_LE.PATH_CONFIG_DEBUG)

    # get the target preparation Enum and prepare a configuration
    try:
        config = TargetPreparationContainer(conf=args.conf, validation=args.validation)
    except Exception as e:
        raise TargetPreparationFailed() from e

    # header: process the header once before actually executing anything
    # ---------
    logger = initialize_logging(config=config, task=_TP.TARGETPREP, _task_enum=_TP, log_conf_path=args.log_conf)
    set_environment(config=config, task=_TP.TARGETPREP, _task_enum=_TP, logger=logger)

    # anything related to Gold (CCDC) fails if the proper environment has not been loaded; now, load the modules here
    try:
        with warnings.catch_warnings(record=True) as w:
            from dockstream.core.Gold.Gold_target_preparator import GoldTargetPreparator
            if len(w) > 0:
                logger.log("Could not load CCDC / Gold target preparator - if another backend is being used, you can safely ignore this warning.", _LE.DEBUG)
    except:
        logger.log("Could not load CCDC / Gold target preparator - if another backend is being used, you can safely ignore this warning.", _LE.DEBUG)

    # make a list of temporary files, that are to be deleted at the end
    temp_files = []

    # do the PDB fixing (if specified)
    # note, that as these steps are backend-independent, we can use the base Enum
    input_pdb_path = config[_TP.TARGETPREP][_TP.INPUT_PATH]
    if config[_TP.TARGETPREP][_TP.FIX][_TP.FIX_ENABLED]:
        pdb_prep = PDBPreparator(conf=config)

        # generate a temporary PDB file, that will be the input later
        temp_pdb_file = gen_temp_file(suffix=".pdb")

        # apply the specified fixing and set the input PDB file
        pdb_prep.fix_pdb(input_pdb_file=input_pdb_path,
                         output_pdb_file=temp_pdb_file)
        temp_files.append(temp_pdb_file)

        # clean-up (and make copy in case specified)
        if nested_get(config, [_TP.TARGETPREP, _TP.FIX, _TP.FIX_PBDOUTPUTPATH], default=False):
            try:
                copyfile(src=temp_pdb_file, dst=config[_TP.TARGETPREP][_TP.FIX][_TP.FIX_PBDOUTPUTPATH])
                input_pdb_path = config[_TP.TARGETPREP][_TP.FIX][_TP.FIX_PBDOUTPUTPATH]
            except:
                logger.log("Could not write fixed intermediate PDB file.", _LE.WARNING)
        else:
            input_pdb_path = temp_pdb_file
        logger.log(f"Wrote fixed PDB to file {input_pdb_path}.", _LE.DEBUG)

    # loop over the specified target preparation steps and execute them
    for run_number, run in enumerate(config[_TP.TARGETPREP][_TP.RUNS]):
        logger.log(f"Started preparation run number {run_number}.", _LE.INFO)
        try:
            if run[_TP.RUNS_BACKEND] == _TP.RUNS_BACKEND_RDOCK:
                prep = rDockTargetPreparator(conf=config, target=input_pdb_path, run_number=run_number)
                result = prep.specify_cavity()
                logger.log("Wrote rDock cavity files to folder specified.",
                           _LE.INFO)
            elif run[_TP.RUNS_BACKEND] == _TP.RUNS_BACKEND_OPENEYE:
                _OpenEye_TP = OpenEyeTargetPreparationEnum()
                prep = OpenEyeTargetPreparator(conf=config, target=input_pdb_path, run_number=run_number)
                prep.specify_cavity()
                prep.write_target(path=run[_TP.RUNS_OUTPUT][_OpenEye_TP.OUTPUT_RECEPTORPATH])
                logger.log(f"Wrote OpenEye receptor to file {run[_TP.RUNS_OUTPUT][_OpenEye_TP.OUTPUT_RECEPTORPATH]}.",
                           _LE.INFO)
            elif run[_TP.RUNS_BACKEND] == _TP.RUNS_BACKEND_GOLD:
                _Gold_TP = GoldTargetPreparationEnum()
                prep = GoldTargetPreparator(conf=config, target=input_pdb_path, run_number=run_number)
                prep.specify_cavity()
                prep.write_target(path=run[_TP.RUNS_OUTPUT][_Gold_TP.OUTPUT_RECEPTORPATH])
                logger.log(f"Wrote Gold target to file {run[_TP.RUNS_OUTPUT][_Gold_TP.OUTPUT_RECEPTORPATH]}.",
                           _LE.INFO)
            elif run[_TP.RUNS_BACKEND] == _TP.RUNS_BACKEND_AUTODOCKVINA:
                _AD_TP = AutodockTargetPreparationEnum()
                prep = AutodockVinaTargetPreparator(conf=config, target=input_pdb_path, run_number=run_number)
                prep.specify_cavity()
                prep.write_target(path=run[_AD_TP.RUNS_OUTPUT][_AD_TP.RECEPTOR_PATH])
                logger.log(f"Wrote AutoDock Vina target to file {run[_AD_TP.RUNS_OUTPUT][_AD_TP.RECEPTOR_PATH]}.",
                           _LE.INFO)
            else:
                raise TargetPreparationFailed("Target preparation backend unknown.")
        except Exception as e:
            logger.log(f"Failed when target preparation run number {run_number}.", _LE.EXCEPTION)
            raise TargetPreparationFailed() from e
        else:
            logger.log(f"Completed target preparation run number {run_number}.", _LE.INFO)

    # clean-up
    for temp_file in temp_files:
        os.remove(temp_file)

    sys.exit(0)
