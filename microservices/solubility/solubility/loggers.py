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
    def starting_api(self):
        _logger.info('Starting api')

    def cuda_available(self, cuda_available):
        _logger.info(f'Cuda available" {str(cuda_available)}')

    def loading_solubility_model(self):
        _logger.info('Loading solubility model')

    def solubility_has_been_loaded(self):
        _logger.info('Solubility model has been loaded')

    def transferred_models_to_gpu(self):
        _logger.info('Solubility model has been loaded')

    def making_solubility_predictions(self):
        _logger.info('Making solubility predictions')

    def successfully_predicted_solubility(self):
        _logger.info('Successfully predicted solubility')

    def run_solubility_prediction_request(self, request: RunSolubilityPredictionRequest):
        d = request.as_log_dict()
        _logger.info('Run solubility request', extra=d)

    def run_solubility_prediction_response(self, response: RunSolubilityPredictionResponse):
        d = response.as_log_dict()
        _logger.info('Run solubility response', extra=d)

    def exception(self):
        _logger.exception('Exception occured in microservice')

logger = Log()