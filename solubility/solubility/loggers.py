import logging
from pythonjsonlogger import jsonlogger

from solubility.api_models import RunSolubilityPredictionRequest, RunSolubilityPredictionResponse

_logger = logging.getLogger()
_logger.setLevel(level=logging.DEBUG)

_logHandler = logging.StreamHandler()
_formatter = jsonlogger.JsonFormatter()
_logHandler.setFormatter(_formatter)
_logger.addHandler(_logHandler)


class Log:
    @staticmethod
    def starting_api():
        _logger.info('Starting api')

    @staticmethod
    def cuda_available(cuda_available):
        _logger.info(f'Cuda available" {str(cuda_available)}')

    @staticmethod
    def loading_solubility_model():
        _logger.info('Loading solubility model')

    @staticmethod
    def solubility_has_been_loaded():
        _logger.info('Solubility model has been loaded')

    @staticmethod
    def transferred_models_to_gpu():
        _logger.info('Solubility model has been loaded')

    @staticmethod
    def making_solubility_predictions():
        _logger.info('Making solubility predictions')

    @staticmethod
    def successfully_predicted_solubility():
        _logger.info('Successfully predicted solubility')

    @staticmethod
    def run_solubility_prediction_request(request: RunSolubilityPredictionRequest):
        d = request.as_log_dict()
        _logger.info('Run solubility request', extra=d)

    @staticmethod
    def run_solubility_prediction_response(response: RunSolubilityPredictionResponse):
        d = response.as_log_dict()
        _logger.info('Run solubility response', extra=d)

    @staticmethod
    def exception():
        _logger.exception('Exception occured in microservice')
