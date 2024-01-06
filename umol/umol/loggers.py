import logging
from pythonjsonlogger import jsonlogger

from umol.api_models import *

_logger = logging.getLogger()

_logHandler = logging.StreamHandler()
_formatter = jsonlogger.JsonFormatter()
_logHandler.setFormatter(_formatter)
_logger.addHandler(_logHandler)


class Log:
    @staticmethod
    def umol_request(request: RunUmolPredictionRequest):
        d = request.as_log_dict()
        _logger.info('Run umol request', extra=d)

    @staticmethod
    def umol_response(response: RunUmolPredictionResponce):
        d = response.as_log_dict()
        _logger.info('Run umol response', extra=d)

    @staticmethod
    def exception():
        _logger.exception()
