from typing import Annotated

from fastapi import FastAPI, File

from microservice.api_models import RunRosettaFoldResponse, IsJobRunningResponse
from microservice.services import RosettaService

app = FastAPI(
    title='RoseTTAFold API',
    version='0.1'
)

from microservice.loggers import logger


@app.post("/run-folding")
async def run_folding(job_id: str | None = None, fasta: Annotated[bytes, File()] = None,
                      a3m: Annotated[bytes, File()] = None) -> RunRosettaFoldResponse:
    """
    Run folding on a given amino-acid sequence
    :param job_id str | None: Used for job tracking
    :param fasta Annotated[bytes, File()]: Fasta file with amino acid sequence. You must specify either fasta or a3m
    :param a3m Annotated[bytes, File()]: MSA file. You must specify either fasta or a3m"""
    if fasta is None and a3m is None:
        return RunRosettaFoldResponse(
            errors=[
                'You must specify either fastsa or a3m'
            ],
            pdb_content=None
        )

    logger.run_rosettafold_request(job_id, fasta, a3m)
    result = await RosettaService().run_rosettafold(job_id, fasta, a3m)
    logger.run_rosettafold_response(result)
    return result


@app.get("/job/{job_id}/is-running")
def is_job_running(job_id: str) -> IsJobRunningResponse:
    return IsJobRunningResponse(is_running=RosettaService().is_job_running(job_id))


@app.get('/livez')
async def livez():
    return True
