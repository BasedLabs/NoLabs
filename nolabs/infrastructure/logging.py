import logging.config
import sys
import traceback

from infrastructure.settings import settings

logger = logging.root

if logger.hasHandlers():
    logger.handlers.clear()

config = {
    'version': 1,
    'disable_existing_loggers': True,
    'loggers': {
        'nolabs': {
            'handlers': ['console'], 'level': settings.logging_level
        }
    },
    'handlers': {
        'console': {
            'formatter': 'nolabs',
            'class': 'logging.StreamHandler'
        }
    },
    'formatters': {
        'nolabs': {
            'format': '%(message)s',
            'class': 'nolabs.infrastructure.logging_formatters.NoLabsFormatter'
        }
    }
}

logging.config.dictConfig(config)
logger = logging.getLogger('nolabs')


def exceptions_hook(exctype, value, tb):
    logger.error(''.join(traceback.format_exception(exctype, value, tb)))


sys.excepthook = exceptions_hook
