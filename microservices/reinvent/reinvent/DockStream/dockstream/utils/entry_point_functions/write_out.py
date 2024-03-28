from dockstream.utils.enums.docking_enum import DockingConfigurationEnum, ResultKeywordsEnum
from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum
from dockstream.utils.enums.logging_enums import LoggingConfigEnum

from dockstream.utils.general_utils import *


def handle_poses_writeout(docking_run, docker, output_prefix):
    _LE = LoggingConfigEnum()
    _LP = LigandPreparationEnum()
    _DE = DockingConfigurationEnum()

    if in_keys(docking_run, [_DE.OUTPUT, _DE.OUTPUT_POSES]):
        if in_keys(docking_run, [_DE.OUTPUT, _DE.OUTPUT_POSES, _DE.OUTPUT_POSES_PATH]):
            poses_path = docking_run[_DE.OUTPUT][_DE.OUTPUT_POSES][_DE.OUTPUT_POSES_PATH]
            poses_path = docker.apply_prefix_to_filename(poses_path, output_prefix)

            # if the overwrite flag is set and the output file exists already, append number to basename
            if nested_get(docking_run, [_DE.OUTPUT,
                                        _DE.OUTPUT_POSES,
                                        _DE.OUTPUT_POSES_OVERWRITE],
                          default=False):
                poses_path = docker.update_path_to_unused(path=poses_path)

            mode = nested_get(docking_run, [_DE.OUTPUT, _DE.OUTPUT_POSES, _DE.OUTPUT_MODE],
                              default=_DE.OUTPUT_MODE_ALL)
            docker.write_docked_ligands(path=poses_path, mode=mode)


def handle_scores_writeout(docking_run, docker, output_prefix):
    _LE = LoggingConfigEnum()
    _LP = LigandPreparationEnum()
    _DE = DockingConfigurationEnum()

    if in_keys(docking_run, [_DE.OUTPUT, _DE.OUTPUT_SCORES]):
        if in_keys(docking_run, [_DE.OUTPUT, _DE.OUTPUT_SCORES, _DE.OUTPUT_SCORES_PATH]):
            scores_path = docking_run[_DE.OUTPUT][_DE.OUTPUT_SCORES][_DE.OUTPUT_SCORES_PATH]
            scores_path = docker.apply_prefix_to_filename(scores_path, output_prefix)

            # if the overwrite flag is set and the output file exists already, append number to basename
            if nested_get(docking_run, [_DE.OUTPUT, _DE.OUTPUT_SCORES, _DE.OUTPUT_SCORES_OVERWRITE],
                          default=False):
                scores_path = docker.update_path_to_unused(path=scores_path)

            mode = nested_get(docking_run, [_DE.OUTPUT, _DE.OUTPUT_SCORES, _DE.OUTPUT_MODE],
                              default=_DE.OUTPUT_MODE_ALL)
            docker.write_result(path=scores_path, mode=mode)


def handle_score_printing(print_scores: bool, print_all: bool, docker, logger):
    _LE = LoggingConfigEnum()
    _LP = LigandPreparationEnum()
    _DE = DockingConfigurationEnum()

    if print_scores:
        _RK = ResultKeywordsEnum()
        scores = docker.get_scores(best_only=not print_all)
        for score in scores:
            print(score, end="\n")
        logger.log(f"Printed {len(scores)} scores to console (print_all set to {print_all}).", _LE.DEBUG)