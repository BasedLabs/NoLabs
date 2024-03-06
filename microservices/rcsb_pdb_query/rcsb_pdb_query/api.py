from fastapi import FastAPI
from rcsb_pdb_query.services import fetch_fasta_files_by_ids, fetch_protein_entries_by_name
from rcsb_pdb_query.api_models import GetFastaFilesByIdsRequest, GetFastaFilesBySearchQueryRequest, \
    GetFastaFilesResponse, IsJobRunningResponse

app = FastAPI(
    title="RCSB PDB Query API"
)

from rcsb_pdb_query.loggers import Log
from rcsb_pdb_query.job_state_manager import job_state_manager

@app.post("/fetch-fastas-by-ids")
def fetch(request: GetFastaFilesByIdsRequest) -> GetFastaFilesResponse:
    if request.job_id:
        job_state_manager.start_job(request.job_id)
    result = fetch_fasta_files_by_ids(request)
    if request.job_id:
        job_state_manager.finish_job(request.job_id)
    return result

@app.post("/fetch-fastas-by-search-query")
def fetch(request: GetFastaFilesBySearchQueryRequest) -> GetFastaFilesResponse:
    if request.job_id:
        job_state_manager.start_job(request.job_id)
    result = fetch_protein_entries_by_name(request)
    if request.job_id:
        job_state_manager.finish_job(request.job_id)
    return result

@app.get("/job/{job_id}/is-running")
def is_job_running(job_id: str) -> IsJobRunningResponse:
    return IsJobRunningResponse(is_running=job_state_manager.is_job_running(job_id))

@app.get("/jobs/running")
def get_running_jobs():
    return {"running_jobs": job_state_manager.get_running_jobs()}
