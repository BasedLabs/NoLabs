import logging
from pythonjsonlogger import jsonlogger

from gene_ontology.api_models import RunGeneOntologyPredictionRequest, RunGeneOntologyPredictionResponse

_logger = logging.getLogger()

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
    def loading_gene_ontology_model():
        _logger.info('Loading gene ontology model')

    @staticmethod
    def gene_ontology_has_been_loaded():
        _logger.info('Gene ontology model has been loaded')

    @staticmethod
    def transferred_models_to_gpu():
        _logger.info('Gene ontology model has been loaded')

    @staticmethod
    def making_gene_ontology_predictions():
        _logger.info('Making gene ontology predictions')

    @staticmethod
    def successfully_predicted_gene_ontology():
        _logger.info('Successfully predicted gene ontology')

    @staticmethod
    def run_go_prediction_request(request: RunGeneOntologyPredictionRequest):
        d = request.as_log_dict()
        _logger.info('Run go request', extra=d)

    @staticmethod
    def run_go_prediction_response(response: RunGeneOntologyPredictionResponse):
        d = response.as_log_dict()
        _logger.info('Run go response', extra=d)

    @staticmethod
    def exception():
        _logger.exception('Exception occured in microservice')
