import pathlib
from typing import Optional, Annotated

from fastapi import FastAPI, UploadFile, File

from microservice.api_models import RunRosettaFoldRequest, RunRosettaFoldResponse
from microservice.services import RosettaService

app = FastAPI(
    title='RoseTTAFold API',
    version='0.1'
)

from microservice.loggers import logger


@app.post("/run-folding")
async def run_folding(fasta: Annotated[bytes, File()] = None, a3m: Annotated[bytes, File()] = None) -> RunRosettaFoldResponse:
    """
    Run folding on a given amino-acid sequence
    :param fasta Annotated[bytes, File()]: Fasta file with amino acid sequence. You must specify either fasta or a3m
    :param a3m Annotated[bytes, File()]: MSA file. You must specify either fasta or a3m"""
    if fasta is None and a3m is None:
        return RunRosettaFoldResponse(
            errors=[
                'You must specify either fastsa or a3m'
            ],
            pdb_content=None
        )

    logger.run_rosettafold_request(fasta, a3m)
    result = await RosettaService().run_rosettafold(fasta, a3m)
    logger.run_rosettafold_response(result)
    return result


@app.get('/instances-running')
async def instances_running() -> int:
    """Get number of rosettafold instances running"""
    return RosettaService().get_rosettafold_process_counter()
