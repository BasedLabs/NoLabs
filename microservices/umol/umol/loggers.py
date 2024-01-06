import logging
from pythonjsonlogger import jsonlogger

from umol.api_models import *

_logger = logging.getLogger()

_logHandler = logging.StreamHandler()
_formatter = jsonlogger.JsonFormatter()
_logHandler.setFormatter(_formatter)
_logger.addHandler(_logHandler)


class Log:
    def umol_request(self, request: RunUmolPredictionRequest):
        d = request.as_log_dict()
        _logger.info('Run umol request', extra=d)

    def umol_response(self, response: RunUmolPredictionResponse):
        d = response.as_log_dict()
        _logger.info('Run umol response', extra=d)

    def exception(self):
        _logger.exception('Exception occured in microservice')


logger = Log()
