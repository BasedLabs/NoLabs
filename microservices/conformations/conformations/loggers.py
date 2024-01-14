import logging
import sys

from pythonjsonlogger import jsonlogger

from conformations.api_models import *

_logger = logging.getLogger()

_logger.setLevel(logging.DEBUG)

_logHandler = logging.StreamHandler(sys.stdout)
_formatter = jsonlogger.JsonFormatter()
_logHandler.setFormatter(_formatter)
_logger.addHandler(_logHandler)


class Log:
    def fixer_request(self, request: RunPdbFixerRequest):
        d = request.as_log_dict()
        _logger.info('Run pdb fixer request', extra=d)

    def fixer_response(self, response: RunPdbFixerResponse):
        d = response.as_log_dict()
        _logger.info('Run pdb fixer response', extra=d)

    def gromacs_simulations_request(self, response: RunGromacsSimulationsRequest):
        d = response.as_log_dict()
        _logger.info('Run gromacs simulations request', extra=d)

    def pdb_simulations_request(self, response: RunPdbSimulationsRequest):
        d = response.as_log_dict()
        _logger.info('Run pdb simulations request', extra=d)

    def gro_top_request(self, response: GenGroTopRequest):
        d = response.as_log_dict()
        _logger.info('Generate gro top request', extra=d)

    def gro_top_response(self, response: GenGroTopResponse):
        d = response.as_log_dict()
        _logger.info('Generate gro top response', extra=d)

    def simulations_response(self, response: RunSimulationsResponse):
        d = response.as_log_dict()
        _logger.info('Run sumulations response', extra=d)

    def exception(self, message: str):
        _logger.exception(message)


logger = Log()
