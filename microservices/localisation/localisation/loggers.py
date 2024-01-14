import logging
from pythonjsonlogger import jsonlogger

from localisation.api_models import RunLocalisationPredictionRequest, RunLocalisationPredictionResponse

_logger = logging.getLogger()
_logger.setLevel(level=logging.DEBUG)

_logHandler = logging.StreamHandler()
_formatter = jsonlogger.JsonFormatter()
_logHandler.setFormatter(_formatter)
_logger.addHandler(_logHandler)


class Log:
    def starting_api(self):
        _logger.info('Starting api')

    def run_localisation_prediction_request(self, request: RunLocalisationPredictionRequest):
        d = request.as_log_dict()
        _logger.info('Run localisation request', extra=d)

    def run_localisation_prediction_response(self, response: RunLocalisationPredictionResponse):
        d = response.as_log_dict()
        _logger.info('Run localisation response', extra=d)


logger = Log()
