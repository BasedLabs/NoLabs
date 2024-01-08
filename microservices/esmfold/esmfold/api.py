from fastapi import FastAPI
from esmfold.services import run_folding, run_facebook_api_folding
from esmfold.api_models import RunEsmFoldPredictionRequest, RunEsmFoldPredictionResponse
from esmfold.loggers import Log

app = FastAPI()

@app.post("/run-folding")
async def predict(request: RunEsmFoldPredictionRequest) -> RunEsmFoldPredictionResponse:
    Log.folding_request(request)
    result = run_folding(request)
    Log.folding_response(result)
    return result

@app.post("/run-through-api-folding")
async def predict_through_api(request: RunEsmFoldPredictionRequest) -> RunEsmFoldPredictionResponse:
    Log.folding_through_api_request(request)
    result = run_facebook_api_folding(request)
    Log.folding_through_api_response(result)
    return result

