import datetime
import json
import logging.handlers
import traceback
from types import TracebackType
from typing import Optional, Tuple, Type, Union

from domain.exceptions import NoLabsException
from nolabs.infrastructure.settings import settings

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

    d = {
        "type": exception_type.__name__,
        "message": str(exception_obj),
        "traceback": (
            "".join(traceback.format_tb(exception_traceback))
            if exception_traceback
            else None
        )
    }

    if exc_info[0] is NoLabsException:
        ex: NoLabsException = exc_info[1]
        if ex.__cause__:
            d['inner_exception'] = ''.join(traceback.TracebackException.from_exception(ex).format())
        d['error_code'] = ex.error_code

    return d


class JsonFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        record_dict = record.__dict__.copy()

        if record_dict.get('args'):
            record_dict['msg'] = record_dict['msg'] % record_dict['args']

            if record_dict.get('colored_message'):
                record_dict['colored_message'] = record_dict['colored_message'] % record_dict['args']

            del record_dict['args']

        record_dict["timestamp"] = (
                datetime.datetime.utcnow().isoformat(timespec="milliseconds") + "Z"
        )
        record_dict["env"] = settings.environment

        record_dict.setdefault("severity", record.levelname)

        if record.exc_info:
            record_dict["exc_info"] = format_exception_info(record.exc_info)

        return json.dumps(record_dict, default=lambda o: str(o))
