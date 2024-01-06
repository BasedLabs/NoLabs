import conformations_microservice as microservice
from conformations_microservice import RunPdbSimulationsRequest, ForceFields, WaterForceFields

import nolabs.api_models.conformations as api
from nolabs.infrastructure.websockets import ConformationsWebsocket
from nolabs.exceptions import NoLabsException, ExceptionCodes


class RunSimulationsFeature:
    def __init__(self, websocket: ConformationsWebsocket):
        self._websocket = websocket

    def _fix_pdb(self, client: microservice.DefaultApi, pdb_content: str, request: api.RunSimulationsRequest) -> str:
        self._websocket.pdb_fixer_start()
        pdb_fixer_request = microservice.RunPdbFixerRequest(
            replace_non_standard_residues=request.replaceNonStandardResidues,
            add_missing_atoms=request.addMissingAtoms,
            add_missing_hydrogens=request.addMissingHydrogens,
            pdb_content=pdb_content
        )
        response = client.fixer_run_pdb_fixer_post(run_pdb_fixer_request=pdb_fixer_request)
        self._websocket.pdb_fixer_finish()
        if response.errors:
            self._websocket.pdb_fixer_errors(response.errors)
            raise NoLabsException(response.errors, error_code=ExceptionCodes.fix_pdb_error)
        return response.pdb_content

    def _run_gromacs_simulations(self, client: microservice.DefaultApi, gro: str, top: str,
                                 request: api.RunSimulationsRequest) -> microservice.RunSimulationsResponse:
        run_gromacs_simulations_request = microservice.RunGromacsSimulationsRequest(
            temperature_k=request.temperatureK,
            friction_coeff=request.frictionCoeff,
            step_size=request.stepSize,
            integrator=microservice.Integrators(request.integrator),
            take_frame_every=request.takeFrameEvery,
            total_frames=request.totalFrames,
            top=top,
            gro=gro
        )
        return client.gromacs_run_gromacs_simulations_post(
            run_gromacs_simulations_request=run_gromacs_simulations_request)

    def _run_pdb_simulations(self, client: microservice.DefaultApi,
                             pdb_content: str,
                             ff: ForceFields,
                             water_ff: WaterForceFields,
                             request: api.RunSimulationsRequest) -> microservice.RunSimulationsResponse:
        run_pdb_simulations_request = RunPdbSimulationsRequest(
            temperature_k=request.temperatureK,
            friction_coeff=request.frictionCoeff,
            step_size=request.stepSize,
            integrator=microservice.Integrators(request.integrator),
            take_frame_every=request.takeFrameEvery,
            total_frames=request.totalFrames,
            pdb_content=pdb_content,
            force_field=ff,
            water_force_field=water_ff
        )
        return client.simulations_run_pdb_simulations_post(run_pdb_simulations_request=run_pdb_simulations_request)

    async def handle(self, request: api.RunSimulationsRequest) -> api.RunSimulationsResponse:
        client = microservice.DefaultApi(api_client=microservice.ApiClient())
        file_content = (await request.pdbFile.read()).decode('utf-8')

        for pdb_content in [file_content, self._fix_pdb(client, file_content, request)]:
            for water_ff in [WaterForceFields(wff) for wff in microservice.WaterForceFields]:
                for ff in [ForceFields(ff) for ff in microservice.ForceFields]:
                    self._websocket.gromacs_gro_top_generator_start(ff, water_ff)
                    generate_gromacs_files_response = client.gen_gro_top_gen_gro_top_post(
                        gen_gro_top_request=microservice.GenGroTopRequest(
                            ignore_missing_atoms=request.ignoreMissingAtoms,
                            force_field=ff,
                            water_force_field=water_ff,
                            pdb_content=pdb_content
                        ))
                    self._websocket.gromacs_simulations_finish()

                    if not generate_gromacs_files_response.errors:
                        self._websocket.gromacs_simulations_start(ff, water_ff)
                        gromacs_simulations_result = self._run_gromacs_simulations(client,
                                                                                   generate_gromacs_files_response.gro,
                                                                                   generate_gromacs_files_response.top,
                                                                                   request)
                        self._websocket.gromacs_simulations_finish()

                        if not gromacs_simulations_result.errors:
                            self._websocket.gromacs_simulations_errors(gromacs_simulations_result.errors)
                            return api.RunSimulationsResponse(gromacs_simulations_result.pdb_content)
                    else:
                        self._websocket.gromacs_gro_top_generator_errors(generate_gromacs_files_response.errors)

                    self._websocket.pdb_simulations_start(ff, water_ff)
                    pdb_simulations_result = self._run_pdb_simulations(client, pdb_content, ff, water_ff, request)
                    self._websocket.pdb_simulations_finish()
                    if not pdb_simulations_result.errors:
                        self._websocket.pdb_simulations_errors(pdb_simulations_result.errors)
                        return api.RunSimulationsResponse(pdb_simulations_result.pdb_content)

        return api.RunSimulationsResponse(errors=['Cannot run simulations for this pdb file'])
