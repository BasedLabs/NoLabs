from fastapi import FastAPI
from solubility.api_models import *
from solubility.services import run_solubility_predictions
from solubility.shared import memory_manager, JobStatusEnum

app = FastAPI(
    title='Solubility api'
)

from solubility.loggers import logger
logger.starting_api()


@app.post("/run")
async def run_solubility(request: RunSolubilityPredictionRequest) -> RunSolubilityPredictionResponse:
    logger.run_solubility_prediction_request(request)
    result = run_solubility_predictions(request)
    logger.run_solubility_prediction_response(result)
    return result


@app.get("/job/{job_id}/is-running")
def is_job_running(job_id: UUID) -> IsJobRunningResponse:
    return IsJobRunningResponse(is_running=memory_manager.get_job_status(job_id) == JobStatusEnum.running)