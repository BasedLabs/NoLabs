import logging
import sys

from pythonjsonlogger import jsonlogger

_logger: logging.Logger = None


def setup_logger(uvicorn_loglevel=logging.ERROR, mongoengine_loglevel=logging.ERROR):
    global _logger
    _logger = logging.getLogger("uvicorn")
    handler = logging.StreamHandler(sys.stdout)
    formatter = jsonlogger.JsonFormatter()

    handler.setFormatter(formatter)
    _logger.addHandler(handler)
    _logger.setLevel(uvicorn_loglevel)

    # Override default handlers for FastAPI
    uvicorn_access_logger = logging.getLogger("uvicorn.access")
    uvicorn_access_logger.addHandler(handler)
    uvicorn_access_logger.setLevel(mongoengine_loglevel)

    uvicorn_error_logger = logging.getLogger("uvicorn.error")
    uvicorn_error_logger.addHandler(handler)
    uvicorn_error_logger.setLevel(mongoengine_loglevel)

    mongoengine_logger = logging.getLogger("mongoengine")
    mongoengine_logger.addHandler(handler)
    mongoengine_logger.setLevel(mongoengine_loglevel)

    mongoengine_logger = logging.getLogger("pymongo.command")
    mongoengine_logger.addHandler(handler)
    mongoengine_logger.setLevel(mongoengine_loglevel)

    mongoengine_logger = logging.getLogger("pymongo")
    mongoengine_logger.addHandler(handler)
    mongoengine_logger.setLevel(mongoengine_loglevel)


def get_logger() -> logging.Logger:
    return _logger
