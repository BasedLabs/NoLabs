from fastapi import FastAPI
app = FastAPI()

from esmfold_light.api_models import RunEsmFoldPredictionRequest, RunEsmFoldPredictionResponse
from esmfold_light.loggers import Log
from esmfold_light.services import run_facebook_api_folding

@app.post("/run-folding")
async def predict_through_api(request: RunEsmFoldPredictionRequest) -> RunEsmFoldPredictionResponse:
    Log.folding_through_api_request(request)
    result = run_facebook_api_folding(request)
    Log.folding_through_api_response(result)
    return result

