from fastapi import FastAPI
from gene_ontology.loggers import logger
from gene_ontology.api_models import *
from gene_ontology.services import run_gene_ontology_prediction

app = FastAPI(
    title='Gene ontology api'
)

logger.starting_api()


@app.post("/run-go-prediction")
async def run_go_prediction(request: RunGeneOntologyPredictionRequest) -> RunGeneOntologyPredictionResponse:
    logger.run_go_prediction_request(request)
    result = run_gene_ontology_prediction(request)
    logger.run_go_prediction_response(result)
    return result
