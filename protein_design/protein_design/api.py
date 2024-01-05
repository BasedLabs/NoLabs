from fastapi import FastAPI
from protein_design.loggers import Log
from protein_design.api_models import *

app = FastAPI()

Log.start_protein_design_api()

@app.post("/run-rfdiffusion")
async def run_rfdiffusion_endpoint(request: RunRfdiffusionRequest) -> RunRfdiffusionResponse:
    Log.run_rfdiffusion_request(request)
    result = run_pdb_fixer(request)
    Log.run_rfdiffusion_response(result)
    return result
