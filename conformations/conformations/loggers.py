import logging
from pythonjsonlogger import jsonlogger

from conformations.api_models import *

_logger = logging.getLogger()

_logHandler = logging.StreamHandler()
_formatter = jsonlogger.JsonFormatter()
_logHandler.setFormatter(_formatter)
_logger.addHandler(_logHandler)


class Log:
    @staticmethod
    def fixer_request(request: RunPdbFixerRequest):
        d = request.as_log_dict()
        _logger.info('Run pdb fixer request', extra=d)

    @staticmethod
    def fixer_response(response: RunPdbFixerResponse):
        d = response.as_log_dict()
        _logger.info('Run pdb fixer response', extra=d)

    @staticmethod
    def gromacs_simulations_request(response: RunGromacsSimulationsRequest):
        d = response.as_log_dict()
        _logger.info('Run gromacs simulations request', extra=d)

    @staticmethod
    def pdb_simulations_request(response: RunPdbSimulationsRequest):
        d = response.as_log_dict()
        _logger.info('Run pdb simulations request', extra=d)

    @staticmethod
    def gro_top_request(response: GenGroTopRequest):
        d = response.as_log_dict()
        _logger.info('Generate gro top request', extra=d)

    @staticmethod
    def gro_top_response(response: GenGroTopResponse):
        d = response.as_log_dict()
        _logger.info('Generate gro top response', extra=d)

    @staticmethod
    def simulations_response(response: RunSimulationsResponse):
        d = response.as_log_dict()
        _logger.info('Run sumulations response', extra=d)

    @staticmethod
    def exception():
        _logger.exception()
