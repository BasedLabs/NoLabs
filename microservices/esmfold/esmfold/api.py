from fastapi import FastAPI
from esmfold.services import run_folding
from esmfold.api_models import RunEsmFoldPredictionRequest, RunEsmFoldPredictionResponse

app = FastAPI(
    title="ESM Fold"
)

from esmfold.loggers import Log


@app.post("/run-folding")
async def predict(request: RunEsmFoldPredictionRequest) -> RunEsmFoldPredictionResponse:
    Log.folding_request(request)
    result = run_folding(request)
    Log.folding_response(result)
    return result

