from conformations.api_models import (GenGroTopRequest, GenGroTopResponse,
                                      IsJobRunningResponse,
                                      RunGromacsSimulationsRequest,
                                      RunPdbFixerRequest, RunPdbFixerResponse,
                                      RunPdbSimulationsRequest,
                                      RunSimulationsResponse)
from conformations.services import (generate_gromacs_files,
                                    run_gromacs_simulation, run_pdb_fixer,
                                    run_pdb_simulation)
from conformations.shared import JobStatusEnum, memory_manager
from fastapi import FastAPI

app = FastAPI(title="Conformations api")

from conformations.loggers import logger


@app.post("/run-pdb-fixer")
async def run_pdb_fixer_endpoint(request: RunPdbFixerRequest) -> RunPdbFixerResponse:
    logger.fixer_request(request)
    result = run_pdb_fixer(request)
    logger.fixer_response(result)
    return result


@app.post("/run-gromacs-simulations")
async def run_gromacs_simulations_endpoint(
    request: RunGromacsSimulationsRequest,
) -> RunSimulationsResponse:
    logger.gromacs_simulations_request(request)
    result = run_gromacs_simulation(request)
    logger.simulations_response(result)
    return result


@app.post("/run-pdb-simulations")
async def run_pdb_simulations_endpoint(
    request: RunPdbSimulationsRequest,
) -> RunSimulationsResponse:
    logger.pdb_simulations_request(request)
    result = run_pdb_simulation(request)
    logger.simulations_response(result)
    return result


@app.post("/gen-gro-top")
async def gen_gro_top_endpoint(request: GenGroTopRequest) -> GenGroTopResponse:
    logger.gro_top_request(request)
    result = generate_gromacs_files(request)
    logger.gro_top_response(result)
    return result


@app.get("/job/{job_id}/is-running")
def is_job_running(job_id: str) -> IsJobRunningResponse:
    return IsJobRunningResponse(
        is_running=memory_manager.get_job_status(job_id) == JobStatusEnum.running
    )
