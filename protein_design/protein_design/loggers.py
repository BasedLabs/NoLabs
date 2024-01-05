import logging
from pythonjsonlogger import jsonlogger

from protein_design.api_models import *

_logger = logging.getLogger()

_logHandler = logging.StreamHandler()
_formatter = jsonlogger.JsonFormatter()
_logHandler.setFormatter(_formatter)
_logger.addHandler(_logHandler)


class Log:
    @staticmethod
    def start_protein_design_api():
        _logger.info('Starting protein design API')

    @staticmethod
    def run_rfdiffusion_request(request: RunRfdiffusionRequest):
        d = request.as_log_dict()
        _logger.info('Run rfdiffusion request', extra=d)

    @staticmethod
    def run_rfdiffusion_response(response: RunRfdiffusionResponse):
        d = response.as_log_dict()
        _logger.info('Run rfdiffusion response', extra=d)

    @staticmethod
    def exception():
        _logger.exception()
