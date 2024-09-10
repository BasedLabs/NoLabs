import logging

from infrastructure.settings import settings

logger = logging.getLogger('nolabs')

logHandler = logging.StreamHandler()
logger.addHandler(logHandler)
logger.setLevel(settings.logging_level)