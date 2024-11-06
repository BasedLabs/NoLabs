import logging
from pythonjsonlogger import jsonlogger

from protmpnn.api_models import *

_logger = logging.getLogger()

_logger.setLevel(logging.DEBUG)

_logHandler = logging.StreamHandler()
_formatter = jsonlogger.JsonFormatter()
_logHandler.setFormatter(_formatter)
_logger.addHandler(_logHandler)


class Log:
    @staticmethod
    def design_request(request: RunProtMPNNPredictionRequest):
        d = request.as_log_dict()
        _logger.info('Run folding on hardware request', extra=d)

    @staticmethod
    def design_response(response: RunProtMPNNPredictionResponse):
        d = response.as_log_dict()
        _logger.info('Run folding on hardware response', extra=d)

    @staticmethod
    def cuda_is_avaialable(availability: bool):
        _logger.info('Cuda is available', extra={
            'availability': availability
        })

    @staticmethod
    def preparing_fasta():
        _logger.info('Preparing fasta file..')

    @staticmethod
    def running_inference():
        _logger.info('Running inference...')


    @staticmethod
    def exception():
        _logger.exception('Exception occured in protmpnn')

logger = Log()