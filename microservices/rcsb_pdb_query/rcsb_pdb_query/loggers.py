import logging
from pythonjsonlogger import jsonlogger

from rcsb_pdb_query.api_models import *

_logger = logging.getLogger()

_logger.setLevel(logging.DEBUG)

_logHandler = logging.StreamHandler()
_formatter = jsonlogger.JsonFormatter()
_logHandler.setFormatter(_formatter)
_logger.addHandler(_logHandler)


class Log:
    @staticmethod
    def query_request(request: GetFastaFilesByIdsRequest):
        d = request.as_log_dict()
        _logger.info('Run request', extra=d)

    @staticmethod
    def query_response(response: GetFastaFilesResponse):
        d = response.as_log_dict()
        _logger.info('Return response', extra=d)

    @staticmethod
    def exception():
        _logger.exception('Exception occured')

logger = Log()