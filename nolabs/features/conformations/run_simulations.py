import conformations_microservice as microservice
from conformations_microservice import RunPdbSimulationsRequest, ForceFields, WaterForceFields

import nolabs.api_models.conformations as api
from nolabs.infrastructure.settings import Settings
from nolabs.features.conformations.services.file_management import FileManagement
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.infrastructure.websockets import ConformationsWebsocket
from nolabs.exceptions import NoLabsException, ErrorCodes


class RunSimulationsFeature:
    def __init__(self, file_management: FileManagement, settings: Settings, websocket: ConformationsWebsocket):
        assert websocket

        self._file_management = file_management
        self._settings = settings
        self._websocket = websocket

    def _fix_pdb(self, client: microservice.DefaultApi, pdb_content: str, request: api.RunSimulationsRequest) -> str:
        self._websocket.pdb_fixer_start()
        pdb_fixer_request = microservice.RunPdbFixerRequest(
            replaceNonStandardResidues=request.replaceNonStandardResidues,
            addMissingAtoms=request.addMissingAtoms,
            addMissingHydrogens=request.addMissingHydrogens,
            pdbContent=pdb_content
        )
        response = client.fixer_run_pdb_fixer_post(run_pdb_fixer_request=pdb_fixer_request)
        self._websocket.pdb_fixer_finish()
        if not response.pdb_content:
            raise NoLabsException('Unexpected error in pdb fixer', error_code=ErrorCodes.fix_pdb_error)
        if response.errors:
            self._websocket.pdb_fixer_errors(response.errors)
            raise NoLabsException(response.errors, error_code=ErrorCodes.fix_pdb_error)
        return response.pdb_content

    def _run_gromacs_simulations(self, client: microservice.DefaultApi, gro: str, top: str,
                                 request: api.RunSimulationsRequest) -> microservice.RunSimulationsResponse:
        run_gromacs_simulations_request = microservice.RunGromacsSimulationsRequest(
            temperatureK=request.temperatureK,
            frictionCoeff=request.frictionCoeff,
            stepSize=request.stepSize,
            integrator=microservice.Integrators(request.integrator),
            takeFrameEvery=request.takeFrameEvery,
            totalFrames=request.totalFrames,
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
            temperatureK=request.temperatureK,
            frictionCoeff=request.frictionCoeff,
            stepSize=request.stepSize,
            integrator=microservice.Integrators(request.integrator),
            takeFrameEvery=request.takeFrameEvery,
            totalFrames=request.totalFrames,
            pdbContent=pdb_content,
            forceField=ff,
            waterForceField=water_ff
        )
        return client.simulations_run_pdb_simulations_post(run_pdb_simulations_request=run_pdb_simulations_request)

    async def handle(self,
                     request: api.RunSimulationsRequest) -> api.RunSimulationsResponse:
        assert request

        experiment_id = ExperimentId(request.experimentId)
        experiment_name = ExperimentName(request.experimentName)

        self._file_management.create_experiment_folder(experiment_id)
        self._file_management.update_metadata(experiment_id, experiment_name, request)

        client = microservice.DefaultApi(api_client=microservice.ApiClient())
        file_content = (await request.pdbFile.read()).decode('utf-8')

        for pdb_content in [file_content, self._fix_pdb(client, file_content, request)]:
            for water_ff in [WaterForceFields(wff) for wff in microservice.WaterForceFields]:
                for ff in [ForceFields(ff) for ff in microservice.ForceFields]:
                    self._websocket.gromacs_gro_top_generator_start(ff, water_ff)
                    generate_gromacs_files_response = client.gen_gro_top_gen_gro_top_post(
                        gen_gro_top_request=microservice.GenGroTopRequest(
                            ignoreMissingAtoms=request.ignoreMissingAtoms,
                            forceField=ff,
                            waterForceField=water_ff,
                            pdbContent=pdb_content
                        ))
                    self._websocket.gromacs_simulations_finish()

                    if not generate_gromacs_files_response.errors and generate_gromacs_files_response.gro and generate_gromacs_files_response.top:
                        self._websocket.gromacs_simulations_start(ff, water_ff)
                        gromacs_simulations_result = self._run_gromacs_simulations(client,
                                                                                   generate_gromacs_files_response.gro,
                                                                                   generate_gromacs_files_response.top,
                                                                                   request)
                        self._websocket.gromacs_simulations_finish()

                        if not gromacs_simulations_result.errors and gromacs_simulations_result.pdb_content:
                            self._websocket.gromacs_simulations_errors(gromacs_simulations_result.errors)

                            self._file_management.save_experiment(experiment_id, gromacs_simulations_result.pdb_content)

                            return api.RunSimulationsResponse(gromacs_simulations_result.pdb_content)
                    else:
                        self._websocket.gromacs_gro_top_generator_errors(generate_gromacs_files_response.errors)

                    self._websocket.pdb_simulations_start(ff, water_ff)
                    pdb_simulations_result = self._run_pdb_simulations(client, pdb_content, ff, water_ff, request)
                    self._websocket.pdb_simulations_finish()
                    if not pdb_simulations_result.errors and pdb_simulations_result.pdb_content:
                        self._websocket.pdb_simulations_errors(pdb_simulations_result.errors)

                        self._file_management.save_experiment(experiment_id, pdb_simulations_result.pdb_content)

                        return api.RunSimulationsResponse(pdb_simulations_result.pdb_content)

        return api.RunSimulationsResponse(errors=['Cannot run simulations for this pdb file'])
