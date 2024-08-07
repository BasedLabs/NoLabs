import os
import ssl
from typing import Dict, Any

from fastapi import FastAPI, HTTPException

from blast_query.job_state_manager import job_state_manager
from blast_query.api_models import SequenceQuery, BlastType, IsJobRunningResponse
from Bio.Blast import NCBIWWW
import xmltodict
import urllib.request

app = FastAPI(title="BLAST Query")


def choose_dataset(program):
    if program in ["blastn", "tblastx", 'tblastn']:
        return "nt"
    elif program in ["blastp", "blastx"]:
        return "nr"
    raise HTTPException(status_code=400, detail="Invalid program")


def run_blast(program: BlastType, query: str, descriptions: int = 10, alignments: int = 10, hitlist_size: int = 10,
              expect: float = 10.0) -> Dict[str, Any]:
    NCBIWWW.email = os.getenv("EMAIL")

    # Create an SSL context that does not verify certificates
    ssl_context = ssl.create_default_context()
    ssl_context.check_hostname = False
    ssl_context.verify_mode = ssl.CERT_NONE

    try:
        dataset = choose_dataset(program)
        # Patch urllib to use the custom SSL context
        opener = urllib.request.build_opener(urllib.request.HTTPSHandler(context=ssl_context))
        urllib.request.install_opener(opener)
        result_handle = NCBIWWW.qblast(program.value, dataset, query, descriptions=descriptions, alignments=alignments,
                                       hitlist_size=hitlist_size, expect=expect)
        result_dict = xmltodict.parse(result_handle.read())
        result_handle.close()
        return result_dict
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/blast")
def blast(query: SequenceQuery) -> Dict[str, Any]:
    """
    Perform a BLAST query with the specified type.

    Args:
        query (SequenceQuery): The query parameters including the sequence and BLAST type.

    Returns:
        Dict[str, Any]: The result of the BLAST query.
    """
    if query.job_id:
        job_state_manager.start_job(query.job_id)
    result = run_blast(query.type, query.sequence, query.descriptions, query.alignments, query.hitlist_size, query.expect)
    if query.job_id:
        job_state_manager.finish_job(query.job_id)
    return result

@app.get("/job/{job_id}/is-running")
def is_job_running(job_id: str) -> IsJobRunningResponse:
    return IsJobRunningResponse(is_running=job_state_manager.is_job_running(job_id))

@app.get("/jobs/running")
def get_running_jobs():
    return {"running_jobs": job_state_manager.get_running_jobs()}