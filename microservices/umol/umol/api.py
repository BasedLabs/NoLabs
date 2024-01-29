from fastapi import FastAPI

from umol.job_state_manager import job_state_manager
from umol.services import run_umol
from umol.api_models import RunUmolPredictionRequest, RunUmolPredictionResponse

app = FastAPI(
    title='Umol'
)

from umol.loggers import logger

@app.post("/run-umol")
async def predict(request: RunUmolPredictionRequest) -> RunUmolPredictionResponse:
    logger.umol_request(request)
    if request.job_id:
        job_state_manager.start_job(request.job_id)
    result = run_umol(request)
    if request.job_id:
        job_state_manager.finish_job(request.job_id)
    logger.umol_response(result)
    return result

