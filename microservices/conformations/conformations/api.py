from fastapi import FastAPI
from conformations.services import run_pdb_fixer, run_gromacs_simulation, run_pdb_simulation, generate_gromacs_files
from conformations.api_models import RunPdbSimulationsRequest, RunGromacsSimulationsRequest, RunSimulationsResponse, \
    RunPdbFixerResponse, RunPdbFixerRequest, GenGroTopRequest, GenGroTopResponse
from conformations.loggers import logger

app = FastAPI(
    title='Conformations api'
)


@app.post("/run-pdb-fixer")
async def run_pdb_fixer_endpoint(request: RunPdbFixerRequest) -> RunPdbFixerResponse:
    logger.fixer_request(request)
    result = run_pdb_fixer(request)
    logger.fixer_response(result)
    return result


@app.post("/run-gromacs-simulations")
async def run_gromacs_simulations_endpoint(request: RunGromacsSimulationsRequest) -> RunSimulationsResponse:
    logger.gromacs_simulations_request(request)
    result = run_gromacs_simulation(request)
    logger.simulations_response(result)
    return result


@app.post('/run-pdb-simulations')
async def run_pdb_simulations_endpoint(request: RunPdbSimulationsRequest) -> RunSimulationsResponse:
    logger.pdb_simulations_request(request)
    result = run_pdb_simulation(request)
    logger.simulations_response(result)
    return result


@app.post('/gen-gro-top')
async def gen_gro_top_endpoint(request: GenGroTopRequest) -> GenGroTopResponse:
    logger.gro_top_request(request)
    result = generate_gromacs_files(request)
    logger.gro_top_response(result)
    return result
