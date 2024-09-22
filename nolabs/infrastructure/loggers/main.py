import logging
import logging.config
import sys
import traceback
import warnings
from pathlib import Path
from typing import Any, Dict

_local_path = Path(__file__).parent


def init_logger(
    file_config=None, level=None, defaults=None, disable_existing_loggers=True
):
    if not file_config:
        file_config = _local_path / "logging.ini"
    logging.config.fileConfig(
        file_config,
        defaults=defaults,
        disable_existing_loggers=disable_existing_loggers,
    )
    if level is None:
        level = "INFO"

    logger = logging.root
    logger.setLevel(level)
    for handler in logger.handlers:
        handler.level = logger.level
    setup_excepthook(logger)


def init_logger_by_dict(config=None) -> None:
    if not config:
        config = make_logging_config()
    logging.config.dictConfig(config)
    setup_excepthook(logging.root)


def setup_excepthook(logger):
    def handler(exctype, value, traceback_):
        logger.error("".join(traceback.format_exception(exctype, value, traceback_)))

    sys.excepthook = handler


def make_logging_config(level=None, disable_existing_loggers=True) -> Dict[str, Any]:
    if level is None:
        level = "INFO"
    return {
        "version": 1,
        "disable_existing_loggers": disable_existing_loggers,
        "root": {"handlers": ["console"], "level": level},
        "handlers": {
            "console": {
                "formatter": "nolabs",
                "class": "logging.StreamHandler",
            }
        },
        "formatters": {
            "nolabs": {
                "format": "%(message)s",
                "class": "nolabs.infrastructure.log_formatters.JsonFormatter",
            }
        },
    }


def get_config(level=None, disable_existing_loggers=True) -> Dict[str, Any]:
    warnings.warn("asdasdasd")
    return make_logging_config(
        level=level, disable_existing_loggers=disable_existing_loggers
    )
