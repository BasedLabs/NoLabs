from fastapi import FastAPI
from chembl_query.services import search_chembl_molecules
from chembl_query.api_models import ChEMBLMoleculeRequest, ChEMBLMoleculeResponse, IsJobRunningResponse

app = FastAPI(
    title="ChemBL Query API"
)

from chembl_query.loggers import Log
from chembl_query.job_state_manager import job_state_manager

@app.post("/query-chembl")
def query(request: ChEMBLMoleculeRequest) -> ChEMBLMoleculeResponse:
    if request.job_id:
        job_state_manager.start_job(request.job_id)
    result = search_chembl_molecules(request)
    if request.job_id:
        job_state_manager.finish_job(request.job_id)
    return result

@app.get("/job/{job_id}/is-running")
def is_job_running(job_id: str) -> IsJobRunningResponse:
    return IsJobRunningResponse(is_running=job_state_manager.is_job_running(job_id))

@app.get("/jobs/running")
def get_running_jobs():
    return {"running_jobs": job_state_manager.get_running_jobs()}
