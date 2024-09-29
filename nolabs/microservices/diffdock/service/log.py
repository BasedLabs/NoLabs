import logging

from celery import signals
from pythonjsonlogger import jsonlogger

handler = logging.StreamHandler()
formatter = jsonlogger.JsonFormatter()
handler.setFormatter(formatter)


@signals.after_setup_logger.connect
def setup_logger(logger, *args, **kwargs):
    """Setup logger to use JSON formatter."""
    logger.handlers.clear()
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)


logger = logging.getLogger("diffdock")
logger.handlers.clear()
logger.addHandler(handler)
logger.setLevel(logging.INFO)
