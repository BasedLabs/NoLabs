import logging
import logging.config
import sys
import traceback


def handler(exctype, value, traceback_):
    logging.root.error("".join(traceback.format_exception(exctype, value, traceback_)))


sys.excepthook = handler


from nolabs.infrastructure.settings import settings

logger = logging.getLogger("nolabs")


def initialize_logging():
    global logger

    if settings.enable_structured_logging:
        from nolabs.infrastructure.log_formatters import JsonFormatter

        try:
            from uvicorn.config import LOGGING_CONFIG

            LOGGING_CONFIG["formatters"]["json"] = {
                "()": JsonFormatter,
            }
        except ModuleNotFoundError:
            pass

        # Update Uvicorn handlers to use the new formatter
        for handler in LOGGING_CONFIG["handlers"].values():
            if "formatter" in handler:
                handler["formatter"] = "json"

        for handler in logging.root.handlers:
            handler.setFormatter(JsonFormatter())

        LOGGING_CONFIG = {
            "version": 1,
            "disable_existing_loggers": False,
            "formatters": {
                "standard": {
                    "format": "%(message)s",
                    "class": "nolabs.infrastructure.log_formatters.JsonFormatter",
                },
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
                },
                "celery": {  # Celery logger
                    "handlers": ["default"],
                    "level": settings.logging_level,
                    "propagate": False,
                },
            },
        }

        logging.config.dictConfig(LOGGING_CONFIG)

    logger = logging.getLogger("nolabs")
