from fastapi import FastAPI

from esmfold_light.job_state_manager import job_state_manager

app = FastAPI()

from esmfold_light.api_models import RunEsmFoldPredictionRequest, RunEsmFoldPredictionResponse, IsJobRunningResponse
from esmfold_light.loggers import Log
from esmfold_light.services import run_facebook_api_folding

@app.post("/run-folding")
def predict_through_api(request: RunEsmFoldPredictionRequest) -> RunEsmFoldPredictionResponse:
    Log.folding_through_api_request(request)
    if request.job_id:
        job_state_manager.start_job(request.job_id)
    result = run_facebook_api_folding(request)
    if request.job_id:
        job_state_manager.finish_job(request.job_id)
    Log.folding_through_api_response(result)
    return result

@app.get("/job/{job_id}/is-running")
def is_job_running(job_id: str) -> IsJobRunningResponse:
    return IsJobRunningResponse(is_running=job_state_manager.is_job_running(job_id))

@app.get("/jobs/running")
def get_running_jobs():
    return {"running_jobs": job_state_manager.get_running_jobs()}