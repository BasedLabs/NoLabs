import logging
import sys

from pythonjsonlogger import jsonlogger

from protein_design.api_models import *

_logger = logging.getLogger()

_logger.setLevel(logging.DEBUG)

_logHandler = logging.StreamHandler(sys.stdout)
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

    def exception(self, exception):
        _logger.exception(exception)


logger = Log()
