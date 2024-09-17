import datetime
import json
import logging.config
import logging.handlers
import sys
import traceback
from types import TracebackType
from typing import Optional, Tuple, Type, Union

from settings import settings

# imported from https://github.com/PrefectHQ/prefect/blob/main/src/prefect/logging/formatters.py

ExceptionInfoType = Union[
    Tuple[Type[BaseException], BaseException, Optional[TracebackType]],
    Tuple[None, None, None],
]


def format_exception_info(exc_info: ExceptionInfoType) -> dict:
    # if sys.exc_info() returned a (None, None, None) tuple,
    # then there's nothing to format
    if exc_info[0] is None:
        return {}

    (exception_type, exception_obj, exception_traceback) = exc_info
    return {
        "type": exception_type.__name__,
        "message": str(exception_obj),
        "traceback": (
            "".join(traceback.format_tb(exception_traceback))
            if exception_traceback
            else None
        ),
    }


# NoLabs structured logging
class JsonFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        record_dict = record.__dict__.copy()

        record_dict["timestamp"] = (
            datetime.datetime.utcnow().isoformat(timespec="milliseconds") + "Z"
        )
        record_dict["env"] = settings.environment

        # GCP severity detection compatibility
        record_dict.setdefault("severity", record.levelname)

        # replace any exception tuples returned by `sys.exc_info()`
        # with a JSON-serializable `dict`.
        if record.exc_info:
            record_dict["exc_info"] = format_exception_info(record.exc_info)

        return json.dumps(record_dict, default=lambda o: str(o))


LOGGING_CONFIG = {
    "version": 1,
    "disable_existing_loggers": True,
    "formatters": {
        "standard": {"format": "%(message)s", "class": "log.JsonFormatter"},
    },
    "handlers": {
        "default": {
            "level": settings.logging_level,
            "formatter": "standard",
            "class": "logging.StreamHandler",
        },
    },
    "loggers": {
        "": {  # root logger
            "handlers": ["default"],
            "level": settings.logging_level,
            "propagate": False,
        }
    },
}

logging.config.dictConfig(LOGGING_CONFIG)

logger = logging.root


def exceptions_hook(exctype, value, tb):
    logger.error("".join(traceback.format_exception(exctype, value, tb)))


sys.excepthook = exceptions_hook
