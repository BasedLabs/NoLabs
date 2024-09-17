from fastapi import FastAPI
from protein_design.api_models import (IsJobRunningResponse,
                                       RunRfdiffusionRequest,
                                       RunRfdiffusionResponse)
from protein_design.services import run_rfdiffusion
from protein_design.shared import JobStatusEnum, memory_manager

app = FastAPI(title="Protein design api")

from protein_design.loggers import logger

logger.start_protein_design_api()


@app.post("/run")
async def run_rfdiffusion_endpoint(
    request: RunRfdiffusionRequest,
) -> RunRfdiffusionResponse:
    logger.run_rfdiffusion_request(request)
    result = run_rfdiffusion(request)
    logger.run_rfdiffusion_response(result)
    return result


@app.get("/job/{job_id}/is-running")
def is_job_running(job_id: str) -> IsJobRunningResponse:
    return IsJobRunningResponse(
        is_running=memory_manager.get_job_status(job_id) == JobStatusEnum.running
    )
