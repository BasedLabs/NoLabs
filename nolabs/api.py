from fastapi import FastAPI

app = FastAPI()


# DO NOT DEFINE ANY MORE GLOBAL VARIABLES HERE!
def register_conformations():
    from nolabs.api_models.conformations import RunSimulationsRequest, RunSimulationsResponse
    from nolabs.features.conformations.run_simulations import RunSimulationsFeature

    @app.post("/run-simulations")
    async def run_simulations(request: RunSimulationsRequest) -> RunSimulationsResponse:
        return RunSimulationsFeature()

    @app.post("/run-gromacs-simulations")
    async def gromacs(request: RunGromacsSimulationsRequest) -> RunSimulationsResponse:
        return run_gromacs_simulation(request)

    @app.post('/run-pdb-simulations')
    async def simulations(request: RunPdbSimulationsRequest) -> RunSimulationsResponse:
        return run_pdb_simulation(request)

    @app.post('/gen-gro-top')
    async def gen_gro_top(request: GenGroTopRequest) -> GenGroTopResponse:
        return generate_gromacs_files(request)
