from fastapi import FastAPI
from solubility.loggers import Log
from solubility.api_models import *
from solubility.services import run_solubility_predictions

app = FastAPI()

Log.starting_api()


@app.post("/run-solubility-prediction")
async def run_solubility(request: RunSolubilityPredictionRequest) -> RunSolubilityPredictionResponse:
    Log.run_solubility_prediction_request(request)
    result = run_solubility_predictions(request)
    Log.run_solubility_prediction_response(result)
    return result
