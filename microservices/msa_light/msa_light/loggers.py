import logging

from msa_light.api_models import (RunMsaPredictionRequest,
                                  RunMsaPredictionResponse)
from pythonjsonlogger import jsonlogger

_logger = logging.getLogger()

_logHandler = logging.StreamHandler()
_formatter = jsonlogger.JsonFormatter()
_logHandler.setFormatter(_formatter)
_logger.addHandler(_logHandler)


class Log:
    @staticmethod
    def msa_request(request: RunMsaPredictionRequest):
        _logger.info("Run msa request")

    @staticmethod
    def msa_response(response: RunMsaPredictionResponse):
        _logger.info("Run msa response")

    @staticmethod
    def exception():
        _logger.exception(msg="Something went wrong")
