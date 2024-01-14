from fastapi import FastAPI
from localisation.api_models import *
from localisation.services import run_localisation

app = FastAPI(
    title='Solubility api'
)

from localisation.loggers import logger
logger.starting_api()


@app.post("/run-localisation-prediction")
async def run_localisation_prediction(request: RunLocalisationPredictionRequest) -> RunLocalisationPredictionResponse:
    logger.run_localisation_prediction_request(request)
    result = run_localisation(request)
    logger.run_localisation_prediction_response(result)
    return result
