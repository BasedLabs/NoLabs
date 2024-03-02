import pathlib

from fastapi import FastAPI, UploadFile

from rosettafold_microservice.api_models import RunRosettaFoldRequest, RunRosettaFoldResponse
from rosettafold_microservice.services import RosettaService

app = FastAPI(
    title='RoseTTAFold API',
    version='0.1'
)

from rosettafold_microservice.loggers import logger


@app.post("/run-folding")
async def run_folding(fasta: UploadFile | None, a3m: UploadFile | None) -> RunRosettaFoldResponse:
    """
    Run folding on a given amino-acid sequence
    :param fasta UploadFile | None: Fasta file with amino acid sequence. You must specify either fasta or a3m
    :param a3m UploadFile | None: MSA file. You must specify either fasta or a3m"""
    if fasta is None and a3m is None:
        return RunRosettaFoldResponse(
            errors=[
                'You must specify either fastsa or a3m'
            ],
            pdb_content=None
        )

    if fasta is not None and pathlib.Path(fasta.filename).suffix not in ['.fasta', '.fa']:
        return RunRosettaFoldResponse(
            errors=[
                'Fasta file must have .fasta or .fa extension'
            ],
            pdb_content=None
        )

    if a3m is not None and pathlib.Path(a3m.filename).suffix not in ['.a3m']:
        return RunRosettaFoldResponse(
            errors=[
                'MSA file must have .a3m extension'
            ],
            pdb_content=None
        )

    await logger.run_rosettafold_request(fasta or a3m)
    result = await RosettaService().run_rosettafold(fasta, a3m)
    logger.run_rosettafold_response(result)
    return result


@app.get('/instances-running')
async def instances_running() -> int:
    """Get number of rosettafold instances running"""
    return RosettaService().get_rosettafold_process_counter()
