from fastapi import FastAPI
from umol.services import run_umol
from umol.api_models import RunUmolPredictionRequest, RunUmolPredictionResponse
from umol.loggers import logger

app = FastAPI(
    title='Umol'
)


@app.post("/run-umol")
async def predict(request: RunUmolPredictionRequest) -> RunUmolPredictionResponse:
    logger.umol_request(request)
    result = run_umol(request)
    logger.umol_response(result)
    return result

