import conformations_microservice as microservice
import nolabs.api_models.conformations as api
from exceptions import NoLabsException, ExceptionCodes


class RunSimulationsFeature:
    def _fix_pdb(self, client: microservice.DefaultApi, pdb_content: str, request: api.RunSimulationsRequest) -> str:
        pdb_fixer_request = microservice.RunPdbFixerRequest(
            replace_non_standard_residues=request.replaceNonStandardResidues,
            add_missing_atoms=request.addMissingAtoms,
            add_missing_hydrogens=request.addMissingHydrogens,
            pdb_content=pdb_content
        )
        response = client.fixer_run_pdb_fixer_post(run_pdb_fixer_request=pdb_fixer_request)
        if response.errors:
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

    async def handle(self, request: api.RunSimulationsRequest) -> api.RunSimulationsResponse:
        client = microservice.DefaultApi(api_client=microservice.ApiClient())
        file_content = (await request.pdbFile.read()).decode('utf-8')

        for pdb_content in [file_content, self._fix_pdb(client, file_content, request)]:
            for water_ff in microservice.WaterForceFields:
                for ff in microservice.ForceFields:
                    generate_gromacs_files_response = client.gen_gro_top_gen_gro_top_post(
                        gen_gro_top_request=microservice.GenGroTopRequest(
                            ignore_missing_atoms=request.ignoreMissingAtoms,
                            force_field=ff,
                            water_force_field=water_ff,
                            pdb_content=pdb_content
                        ))
                    if not generate_gromacs_files_response.errors:
                        gromacs_simulations_result = self._run_gromacs_simulations(client,
                                                                                   generate_gromacs_files_response.gro,
                                                                                   generate_gromacs_files_response.top,
                                                                                   request)
                        if not gromacs_simulations_result.errors:
                            return api.RunSimulationsResponse(gromacs_simulations_result.pdb_content)

                    client.


        response = client.gromacs_run_gromacs_simulations_post(
            run_gromacs_simulations_request=Run)
        return response
