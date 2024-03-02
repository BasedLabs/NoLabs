from fastapi import FastAPI
from pubmed_query.services import search_pubmed
from pubmed_query.api_models import PubMedSearchRequest, PubMedSearchResponse, IsJobRunningResponse

app = FastAPI(
    title="Pubmed Query API"
)

from pubmed_query.loggers import (Log)
from pubmed_query.job_state_manager import job_state_manager


@app.post("/search_pubmed_articles")
def search(request: PubMedSearchRequest) -> PubMedSearchResponse:
    if request.job_id:
        job_state_manager.start_job(request.job_id)
    result = search_pubmed(request)
    if request.job_id:
        job_state_manager.finish_job(request.job_id)
    return result


@app.get("/job/{job_id}/is-running")
def is_job_running(job_id: str) -> IsJobRunningResponse:
    return IsJobRunningResponse(is_running=job_state_manager.is_job_running(job_id))


@app.get("/jobs/running")
def get_running_jobs():
    return {"running_jobs": job_state_manager.get_running_jobs()}
