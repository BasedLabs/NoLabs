from fastapi import FastAPI
from localisation.api_models import (IsJobRunningResponse,
                                     RunLocalisationPredictionRequest,
                                     RunLocalisationPredictionResponse)
from localisation.services import run_localisation
from localisation.shared import JobStatusEnum, memory_manager

app = FastAPI(title="Solubility api")

from localisation.loggers import logger

logger.starting_api()


@app.post("/run")
async def predict(
    request: RunLocalisationPredictionRequest,
) -> RunLocalisationPredictionResponse:
    logger.run_localisation_prediction_request(request)
    result = run_localisation(request)
    logger.run_localisation_prediction_response(result)
    return result


@app.get("/job/{job_id}/is-running")
def is_job_running(job_id: str) -> IsJobRunningResponse:
    return IsJobRunningResponse(
        is_running=memory_manager.get_job_status(job_id) == JobStatusEnum.running
    )
