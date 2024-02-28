from fastapi import FastAPI
from biobuddy.services import send_message
from biobuddy.api_models import SendMessageToBioBuddyRequest, SendMessageToBioBuddyResponse, IsJobRunningResponse

app = FastAPI(
    title="Bio Buddy"
)

from biobuddy.loggers import Log
from biobuddy.job_state_manager import job_state_manager

@app.post("/send-message")
def predict(request: SendMessageToBioBuddyRequest) -> SendMessageToBioBuddyResponse:
    if request.job_id:
        job_state_manager.start_job(request.job_id)
    result = send_message(request)
    if request.job_id:
        job_state_manager.finish_job(request.job_id)
    return result

@app.get("/job/{job_id}/is-running")
def is_job_running(job_id: str) -> IsJobRunningResponse:
    return IsJobRunningResponse(is_running=job_state_manager.is_job_running(job_id))

@app.get("/jobs/running")
def get_running_jobs():
    return {"running_jobs": job_state_manager.get_running_jobs()}
