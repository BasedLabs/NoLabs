from fastapi import FastAPI
from conformations.services import run_pdb_fixer, run_gromacs_simulation, run_pdb_simulation, generate_gromacs_files
from conformations.api_models import RunPdbSimulationsRequest, RunGromacsSimulationsRequest, RunSimulationsResponse, \
    RunPdbFixerResponse, RunPdbFixerRequest, GenGroTopRequest, GenGroTopResponse
from conformations.loggers import Log

app = FastAPI()


@app.post("/run-pdb-fixer")
async def run_pdb_fixer_endpoint(request: RunPdbFixerRequest) -> RunPdbFixerResponse:
    Log.fixer_request(request)
    result = run_pdb_fixer(request)
    Log.fixer_response(result)
    return result


@app.post("/run-gromacs-simulations")
async def run_gromacs_simulations_endpoint(request: RunGromacsSimulationsRequest) -> RunSimulationsResponse:
    Log.gromacs_simulations_request(request)
    result = run_gromacs_simulation(request)
    Log.simulations_response(result)
    return result


@app.post('/run-pdb-simulations')
async def run_pdb_simulations_endpoint(request: RunPdbSimulationsRequest) -> RunSimulationsResponse:
    Log.pdb_simulations_request(request)
    result = run_pdb_simulation(request)
    Log.simulations_response(result)
    return result


@app.post('/gen-gro-top')
async def gen_gro_top_endpoint(request: GenGroTopRequest) -> GenGroTopResponse:
    Log.gro_top_request(request)
    result = generate_gromacs_files(request)
    Log.gro_top_response(result)
    return result
