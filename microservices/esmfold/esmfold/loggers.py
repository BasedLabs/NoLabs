import logging
from pythonjsonlogger import jsonlogger

from esmfold.api_models import *

_logger = logging.getLogger()

_logHandler = logging.StreamHandler()
_formatter = jsonlogger.JsonFormatter()
_logHandler.setFormatter(_formatter)
_logger.addHandler(_logHandler)


class Log:
    @staticmethod
    def folding_request(request: RunEsmFoldPredictionRequest):
        d = request.as_log_dict()
        _logger.info('Run folding on hardware request', extra=d)

    @staticmethod
    def folding_response(response: RunEsmFoldPredictionResponce):
        d = response.as_log_dict()
        _logger.info('Run folding on hardware response', extra=d)

    @staticmethod
    def folding_through_api_request(request: RunEsmFoldPredictionRequest):
        d = request.as_log_dict()
        _logger.info('Run folding via facebook api request', extra=d)

    @staticmethod
    def folding_through_api_response(response: RunEsmFoldPredictionResponce):
        d = response.as_log_dict()
        _logger.info('Run folding via facebook api response', extra=d)

    @staticmethod
    def exception():
        _logger.exception()
