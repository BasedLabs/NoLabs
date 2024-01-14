import logging
from pythonjsonlogger import jsonlogger

from protein_design.api_models import *

_logger = logging.getLogger()

_logger.setLevel(logging.DEBUG)

_logHandler = logging.StreamHandler()
_formatter = jsonlogger.JsonFormatter()
_logHandler.setFormatter(_formatter)
_logger.addHandler(_logHandler)


class Log:
    def start_protein_design_api(self):
        _logger.info('Starting protein design API')

    def run_rfdiffusion_request(self, request: RunRfdiffusionRequest):
        d = request.as_log_dict()
        _logger.info('Run rfdiffusion request', extra=d)

    def run_rfdiffusion_response(self, response: RunRfdiffusionResponse):
        d = response.as_log_dict()
        _logger.info('Run rfdiffusion response', extra=d)

    def exception(self):
        _logger.exception('Exception occured in microservice')


logger = Log()
