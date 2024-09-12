import sys
import traceback

import logging.config

from infrastructure.settings import settings

LOGGING_CONFIG = {
    'version': 1,
    'disable_existing_loggers': True,
    'formatters': {
        'standard': {
            'format': '%(message)s',
            'class': 'nolabs.infrastructure.log_formatters.JsonFormatter'
        },
    },
    'handlers': {
        'default': {
            'level': settings.logging_level,
            'formatter': 'standard',
            'class': 'logging.StreamHandler'
        },
    },
    'loggers': {
        '': {  # root logger
            'handlers': ['default'],
            'level': settings.logging_level,
            'propagate': False
        }
    }
}

logging.config.dictConfig(LOGGING_CONFIG)

logger = logging.getLogger("nolabs")


def exceptions_hook(exctype, value, tb):
    logger.error(''.join(traceback.format_exception(exctype, value, tb)))


sys.excepthook = exceptions_hook
