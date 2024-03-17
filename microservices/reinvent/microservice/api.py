import pathlib
from typing import Optional, Annotated

from fastapi import FastAPI, File, APIRouter

from microservice.api_models import RunRosettaFoldResponse, IsJobRunningResponse, RunFunTuningJobRequest, \
    RunFineTuningJobRequest, FineTuningJobResponse
from microservice.services import RosettaService

app = FastAPI(
    title='RoseTTAFold API',
    version='0.1'
)

from microservice.loggers import logger


router = APIRouter(
    prefix='/api/v1/fine-tuning',
    tags=['folding']
)


@app.post("/fine-tuning/run")
def run_fine_tuning(request: RunFineTuningJobRequest, pdb_content: Annotated[bytes, File()]) -> FineTuningJobResponse:


