import logging
from pythonjsonlogger import jsonlogger

from msa_light.api_models import *

_logger = logging.getLogger()

_logHandler = logging.StreamHandler()
_formatter = jsonlogger.JsonFormatter()
_logHandler.setFormatter(_formatter)
_logger.addHandler(_logHandler)


class Log:
    @staticmethod
    def msa_request(request: RunMsaPredictionRequest):
        d = request
        _logger.info('Run msa request')

    @staticmethod
    def msa_response(response: RunMsaPredictionResponse):
        d = response
        _logger.info('Run msa response')

    @staticmethod
    def exception():
        _logger.exception(msg='Something went wrong')
