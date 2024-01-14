from fastapi import FastAPI
from solubility.api_models import *
from solubility.services import run_solubility_predictions

app = FastAPI(
    title='Solubility api'
)

from solubility.loggers import logger
logger.starting_api()


@app.post("/run-solubility-prediction")
async def run_solubility(request: RunSolubilityPredictionRequest) -> RunSolubilityPredictionResponse:
    logger.run_solubility_prediction_request(request)
    result = run_solubility_predictions(request)
    logger.run_solubility_prediction_response(result)
    return result
