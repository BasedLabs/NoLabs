from typing import Dict, List

import conformations_microservice as microservice
from conformations_microservice import RunPdbSimulationsRequest, OpenMmForceFields, OpenMmWaterForceFields, \
    GromacsWaterForceFields, \
    GromacsForceFields

import nolabs.api_models.conformations as api
from nolabs.features.events_queue import EventsQueue, EventsQueueMessageClass
from nolabs.infrastructure.settings import Settings
from nolabs.features.conformations.services.file_management import FileManagement
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.utils import generate_uuid


class RunSimulationsFeature:
    def __init__(self,
                 file_management: FileManagement,
                 settings: Settings,
                 events_queue: EventsQueue):
        self._file_management = file_management
        self._settings = settings
        self._events_queue = events_queue

    async def handle(self,
                     request: api.RunSimulationsRequest) -> api.RunSimulationsResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id) if request.experiment_id else ExperimentId(generate_uuid())
        experiment_name = ExperimentName(request.experiment_name)

        self._file_management.create_experiment_folder(experiment_id)
        self._file_management.update_metadata(experiment_id, experiment_name, request)

        configuration = microservice.Configuration(
            host=self._settings.conformations_host,
        )
        with microservice.ApiClient(configuration=configuration) as _:
            client = microservice.DefaultApi(api_client=_)
            file_content = (await request.pdb_file.read()).decode('utf-8')

            for pdb_content in [file_content, self._fix_pdb(client, file_content, request)]:
                for water_ff in [GromacsWaterForceFields(wff) for wff in microservice.GromacsWaterForceFields]:
                    for ff in [GromacsForceFields(ff) for ff in microservice.GromacsForceFields]:
                        self._enqueue_event(
                            f'Gromacs gro top generator start, force field: {ff}, water force field: {water_ff}')
                        generate_gromacs_files_response = client.gen_gro_top_endpoint_gen_gro_top_post(
                            gen_gro_top_request=microservice.GenGroTopRequest(
                                ignore_missing_atoms=request.ignore_missing_atoms,
                                force_field=ff,
                                water_force_field=water_ff,
                                pdb_content=pdb_content
                            ))

                        if generate_gromacs_files_response.errors:
                            self._enqueue_event(message='Gromacs gro top generator errors',
                                                errors=generate_gromacs_files_response.errors)

                        if generate_gromacs_files_response.gro and generate_gromacs_files_response.top:
                            f'Gromacs simulations started, force field: {ff}, water force field: {water_ff}'
                            gromacs_simulations_result = self._run_gromacs_simulations(client,
                                                                                       generate_gromacs_files_response.gro,
                                                                                       generate_gromacs_files_response.top,
                                                                                       request)

                            if gromacs_simulations_result.errors:
                                self._enqueue_event(message='Gromacs simulations errors',
                                                    errors=gromacs_simulations_result.errors)

                            if gromacs_simulations_result.pdb_content:
                                self._file_management.save_experiment(experiment_id,
                                                                      gromacs_simulations_result.pdb_content)
                                self._enqueue_event('Gromacs simulations finished')
                                return api.RunSimulationsResponse(experiment_id=experiment_id.value,
                                                                  pdb_content=gromacs_simulations_result.pdb_content)

                for open_mm_water_ff in [OpenMmWaterForceFields(wff) for wff in microservice.OpenMmWaterForceFields]:
                    for open_mm_ff in [OpenMmForceFields(ff) for ff in microservice.OpenMmForceFields]:
                        self._enqueue_event(
                            message=f'Pdb simulations started, force field: {open_mm_ff}, water force field: {open_mm_water_ff}')
                        pdb_simulations_result = self._run_pdb_simulations(client, pdb_content, open_mm_ff,
                                                                           open_mm_water_ff, request)
                        self._enqueue_event('Pdb simulations finished')

                        if pdb_simulations_result.errors:
                            self._enqueue_event(message='Pdb simulations errors', errors=pdb_simulations_result.errors)

                        if pdb_simulations_result.pdb_content:
                            self._file_management.save_experiment(experiment_id, pdb_simulations_result.pdb_content)

                            return api.RunSimulationsResponse(experiment_id=experiment_id.value,
                                                              pdb_content=pdb_simulations_result.pdb_content)

            return api.RunSimulationsResponse(experiment_id=experiment_id.value,
                                              errors=['Cannot run simulations for this pdb file'])

    def _enqueue_event(self, message: str | None = None, errors: List[str] | None = None):
        j = {
            'message': message,
            'errors': errors
        }
        self._events_queue.send_json(EventsQueueMessageClass.conformations, j)

    def _fix_pdb(self, client: microservice.DefaultApi, pdb_content: str, request: api.RunSimulationsRequest) -> str:
        self._enqueue_event(message='Pdb fixer started')
        pdb_fixer_request = microservice.RunPdbFixerRequest(
            replace_nonstandard_residues=request.replace_non_standard_residues,
            add_missing_atoms=request.add_missing_atoms,
            add_missing_hydrogens=request.add_missing_hydrogens,
            pdb_content=pdb_content
        )
        response = client.run_pdb_fixer_endpoint_run_pdb_fixer_post(run_pdb_fixer_request=pdb_fixer_request)
        self._enqueue_event(message='Pdb fixer finished')
        if not response.pdb_content:
            raise NoLabsException('Unexpected error in pdb fixer', error_code=ErrorCodes.fix_pdb_error)
        if response.errors:
            self._enqueue_event(message='Pdb fixer errors', errors=response.errors)
            raise NoLabsException(response.errors, error_code=ErrorCodes.fix_pdb_error)
        return response.pdb_content

    def _run_gromacs_simulations(self, client: microservice.DefaultApi, gro: str, top: str,
                                 request: api.RunSimulationsRequest) -> microservice.RunSimulationsResponse:
        run_gromacs_simulations_request = microservice.RunGromacsSimulationsRequest(
            temperature_k=request.temperature_k,
            friction_coeff=request.friction_coeff,
            step_size=request.step_size,
            integrator=microservice.Integrators(request.integrator.value),  # type: ignore
            take_frame_every=request.take_frame_every,
            total_frames=request.total_frames,
            top=top,
            gro=gro
        )
        return client.run_gromacs_simulations_endpoint_run_gromacs_simulations_post(
            run_gromacs_simulations_request=run_gromacs_simulations_request)

    def _run_pdb_simulations(self, client: microservice.DefaultApi,
                             pdb_content: str,
                             ff: OpenMmForceFields,
                             water_ff: OpenMmWaterForceFields,
                             request: api.RunSimulationsRequest) -> microservice.RunSimulationsResponse:
        run_pdb_simulations_request = RunPdbSimulationsRequest(
            temperature_k=request.temperature_k,
            friction_coeff=request.friction_coeff,
            step_size=request.step_size,
            integrator=microservice.Integrators(request.integrator.value),  # type: ignore
            take_frame_every=request.take_frame_every,
            total_frames=request.total_frames,
            pdb_content=pdb_content,
            force_field=ff,
            water_force_field=water_ff
        )
        return client.run_pdb_simulations_endpoint_run_pdb_simulations_post(
            run_pdb_simulations_request=run_pdb_simulations_request)
