import logging
from pythonjsonlogger import jsonlogger

from gene_ontology.api_models import RunGeneOntologyPredictionRequest, RunGeneOntologyPredictionResponse

_logger = logging.getLogger()

_logger.setLevel(logging.DEBUG)

_logHandler = logging.StreamHandler()
_formatter = jsonlogger.JsonFormatter()
_logHandler.setFormatter(_formatter)
_logger.addHandler(_logHandler)


class Log:
    def starting_api(self):
        _logger.info('Starting api')

    def cuda_available(self, cuda_available):
        _logger.info(f'Cuda available" {str(cuda_available)}')

    def loading_gene_ontology_model(self):
        _logger.info('Loading gene ontology model')

    def gene_ontology_has_been_loaded(self):
        _logger.info('Gene ontology model has been loaded')

    def transferred_models_to_gpu(self):
        _logger.info('Gene ontology model has been loaded')

    def making_gene_ontology_predictions(self):
        _logger.info('Making gene ontology predictions')

    def successfully_predicted_gene_ontology(self):
        _logger.info('Successfully predicted gene ontology')

    def run_go_prediction_request(self, request: RunGeneOntologyPredictionRequest):
        d = request.as_log_dict()
        _logger.info('Run go request', extra=d)

    def run_go_prediction_response(self, response: RunGeneOntologyPredictionResponse):
        d = response.as_log_dict()
        _logger.info('Run go response', extra=d)

    def exception(self):
        _logger.exception('Exception occured in microservice')

logger = Log()