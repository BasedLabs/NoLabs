from fastapi import FastAPI
from msa_light.job_state_manager import job_state_manager

app = FastAPI()

from msa_light.api_models import (IsJobRunningResponse,
                                  RunMsaPredictionRequest,
                                  RunMsaPredictionResponse)
from msa_light.loggers import Log
from msa_light.services import predict_msa_service


@app.post("/predict-msa")
def predict_msa(request: RunMsaPredictionRequest) -> RunMsaPredictionResponse:
    Log.msa_request(request)
    if request.job_id:
        job_state_manager.start_job(request.job_id)
    result = predict_msa_service(request)
    if request.job_id:
        job_state_manager.finish_job(request.job_id)
    Log.msa_response(result)
    return result


@app.get("/job/{job_id}/is-running")
def is_job_running(job_id: str) -> IsJobRunningResponse:
    return IsJobRunningResponse(is_running=job_state_manager.is_job_running(job_id))


@app.get("/jobs/running")
def get_running_jobs():
    return {"running_jobs": job_state_manager.get_running_jobs()}
