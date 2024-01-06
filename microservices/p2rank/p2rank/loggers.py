import logging
from pythonjsonlogger import jsonlogger

from p2rank.api_models import *

_logger = logging.getLogger()

_logHandler = logging.StreamHandler()
_formatter = jsonlogger.JsonFormatter()
_logHandler.setFormatter(_formatter)
_logger.addHandler(_logHandler)


class Log:
    def p2rank_request(self, request: RunP2RankPredictionRequest):
        d = request.as_log_dict()
        _logger.info('Run p2rank request', extra=d)

    def p2rank_response(self, response: RunP2RankPredictionResponse):
        d = response.as_log_dict()
        _logger.info('Run p2rank response', extra=d)

    def exception(self):
        _logger.exception('Exception occured in microservice')


logger = Log()
