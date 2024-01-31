from fastapi import FastAPI

from p2rank.job_state_manager import job_state_manager
from p2rank.services import run_p2rank
from p2rank.api_models import RunP2RankPredictionRequest, RunP2RankPredictionResponse, IsJobRunningResponse

app = FastAPI()

from p2rank.loggers import logger

@app.post("/run-p2rank")
def predict(request: RunP2RankPredictionRequest) -> RunP2RankPredictionResponse:
    logger.p2rank_request(request)
    if request.job_id:
        job_state_manager.start_job(request.job_id)
    result = run_p2rank(request)
    if request.job_id:
        job_state_manager.finish_job(request.job_id)
    logger.p2rank_response(result)
    return result

@app.get("/job/{job_id}/is-running")
def is_job_running(job_id: str) -> IsJobRunningResponse:
    return IsJobRunningResponse(is_running=job_state_manager.is_job_running(job_id))

@app.get("/jobs/running")
def get_running_jobs():
    return {"running_jobs": job_state_manager.get_running_jobs()}

