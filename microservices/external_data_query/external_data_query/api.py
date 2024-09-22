from typing import Any, Dict

from external_data_query.chembl.api_models import (ChEMBLMoleculeRequest,
                                                   ChEMBLMoleculeResponse,
                                                   DrugIndicationRequest,
                                                   DrugIndicationResponse)
from external_data_query.chembl.services import (search_chembl_molecules,
                                                 search_drugs_for_condition)
from external_data_query.pubmed.api_models import (PubMedSearchRequest,
                                                   PubMedSearchResponse)
from external_data_query.pubmed.services import search_pubmed
from external_data_query.rcsb_pdb.api_models import (
    ComplexQueryRequest, GetFastaFilesByIdsRequest,
    GetFastaFilesBySearchQueryRequest, GetFastaFilesResponse,
    IsJobRunningResponse, SequenceQueryRequest)
from external_data_query.rcsb_pdb.services import (
    create_query_node, fetch_entries_by_complex_query,
    fetch_fasta_files_by_ids, fetch_protein_entries_by_name,
    fetch_proteins_by_sequence)
from fastapi import FastAPI

app = FastAPI(title="External Databases Query API")

from external_data_query.job_state_manager import job_state_manager


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


@app.post("/fetch-fastas-by-complex-query")
def fetch(request: Dict[str, Any]) -> GetFastaFilesResponse:
    request = ComplexQueryRequest(
        query=create_query_node(request["query"]),
        max_results=request.get("max_results"),
        job_id=request.get("job_id"),
    )
    if request.job_id:
        job_state_manager.start_job(request.job_id)
    result = fetch_entries_by_complex_query(request)
    if request.job_id:
        job_state_manager.finish_job(request.job_id)
    return result


@app.post("/fetch-fastas-by-sequence")
def fetch(request: SequenceQueryRequest) -> GetFastaFilesResponse:
    if request.job_id:
        job_state_manager.start_job(request.job_id)
    result = fetch_proteins_by_sequence(request)
    if request.job_id:
        job_state_manager.finish_job(request.job_id)
    return result


@app.post("/query-chembl")
def query(request: ChEMBLMoleculeRequest) -> ChEMBLMoleculeResponse:
    if request.job_id:
        job_state_manager.start_job(request.job_id)
    result = search_chembl_molecules(request)
    if request.job_id:
        job_state_manager.finish_job(request.job_id)
    return result


@app.post("/query-chembl-by-condition")
def query(request: DrugIndicationRequest) -> DrugIndicationResponse:
    if request.job_id:
        job_state_manager.start_job(request.job_id)
    result = search_drugs_for_condition(request)
    if request.job_id:
        job_state_manager.finish_job(request.job_id)
    return result


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
