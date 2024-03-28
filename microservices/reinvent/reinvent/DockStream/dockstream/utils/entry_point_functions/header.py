import os

import logging.config as logging_config
from dockstream.loggers.interface_logger import InterfaceLogger

from dockstream.utils.files_paths import dict_from_json_file

from dockstream.utils.enums.logging_enums import LoggingConfigEnum
from dockstream.utils.general_utils import *

_LE = LoggingConfigEnum()


def initialize_logging(config, task, _task_enum, log_conf_path):
    if in_keys(config, [task, _task_enum.HEADER]):
        log_conf_dict = dict_from_json_file(log_conf_path)
        if in_keys(config, [task, _task_enum.HEADER, _task_enum.LOGGING]):
            if in_keys(config, [task, _task_enum.HEADER, _task_enum.LOGGING, _task_enum.LOGGING_LOGFILE]):
                try:
                    log_conf_dict["handlers"]["file_handler"]["filename"] = config[task][_task_enum.HEADER][_task_enum.LOGGING][_task_enum.LOGGING_LOGFILE]
                    log_conf_dict["handlers"]["file_handler_blank"]["filename"] = config[task][_task_enum.HEADER][_task_enum.LOGGING][_task_enum.LOGGING_LOGFILE]
                except KeyError:
                    pass
        logging_config.dictConfig(log_conf_dict)
    else:
        logging_config.dictConfig(dict_from_json_file(log_conf_path))
    logger = InterfaceLogger()
    logger.log(f"DockStream version used: {parse_setuppy()['version']}", _LE.INFO)
    return logger


def set_environment(config, task, _task_enum, logger):
    if in_keys(config, [task, _task_enum.HEADER]):
        if in_keys(config, [task, _task_enum.HEADER, _task_enum.ENVIRONMENT]):
            if in_keys(config, [task, _task_enum.HEADER, _task_enum.ENVIRONMENT, _task_enum.ENVIRONMENT_EXPORT]):
                exp_vars = config[task][_task_enum.HEADER][_task_enum.ENVIRONMENT][_task_enum.ENVIRONMENT_EXPORT]
                for export in exp_vars:
                    os.environ[export[_task_enum.ENVIRONMENT_EXPORT_KEY]] = export[_task_enum.ENVIRONMENT_EXPORT_VALUE]
                    logger.log(f"Added environment variable {export[_task_enum.ENVIRONMENT_EXPORT_KEY]}: {export[_task_enum.ENVIRONMENT_EXPORT_VALUE]}.", _LE.DEBUG)
