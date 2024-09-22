from esmfold.api_models import (IsJobRunningResponse,
                                RunEsmFoldPredictionRequest,
                                RunEsmFoldPredictionResponse)
from esmfold.services import run_folding
from fastapi import FastAPI

app = FastAPI(title="ESM Fold")

from esmfold.job_state_manager import job_state_manager
from esmfold.loggers import Log


@app.post("/run-folding")
def predict(request: RunEsmFoldPredictionRequest) -> RunEsmFoldPredictionResponse:
    Log.folding_request(request)
    if request.job_id:
        job_state_manager.start_job(request.job_id)
    result = run_folding(request)
    if request.job_id:
        job_state_manager.finish_job(request.job_id)
    Log.folding_response(result)
    return result


@app.get("/job/{job_id}/is-running")
def is_job_running(job_id: str) -> IsJobRunningResponse:
    return IsJobRunningResponse(is_running=job_state_manager.is_job_running(job_id))


@app.get("/jobs/running")
def get_running_jobs():
    return {"running_jobs": job_state_manager.get_running_jobs()}
