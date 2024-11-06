from fastapi import FastAPI
from protmpnn.services import run_protmpnn
from protmpnn.api_models import RunProtMPNNPredictionRequest, RunProtMPNNPredictionResponse, IsJobRunningResponse

app = FastAPI(
    title="ProtMPNN"
)

from protmpnn.loggers import Log
from protmpnn.job_state_manager import job_state_manager

@app.post("/run-design")
def predict(request: RunProtMPNNPredictionRequest) -> RunProtMPNNPredictionResponse:
    if request.job_id:
        job_state_manager.start_job(request.job_id)
    result = run_protmpnn(request)
    if request.job_id:
        job_state_manager.finish_job(request.job_id)
    return result

@app.get("/job/{job_id}/is-running")
def is_job_running(job_id: str) -> IsJobRunningResponse:
    return IsJobRunningResponse(is_running=job_state_manager.is_job_running(job_id))

@app.get("/jobs/running")
def get_running_jobs():
    return {"running_jobs": job_state_manager.get_running_jobs()}
