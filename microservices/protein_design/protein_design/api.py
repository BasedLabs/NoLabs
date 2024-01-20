from fastapi import FastAPI
from protein_design.api_models import *

app = FastAPI(
    title='Protein design api'
)

from protein_design.loggers import logger
logger.start_protein_design_api()


@app.post("/run-rfdiffusion")
async def run_rfdiffusion_endpoint(request: RunRfdiffusionRequest) -> RunRfdiffusionResponse:
    logger.run_rfdiffusion_request(request)
    result = await run_rfdiffusion_endpoint(request)
    logger.run_rfdiffusion_response(result)
    return result
