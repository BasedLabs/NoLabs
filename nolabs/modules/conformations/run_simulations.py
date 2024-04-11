from typing import List

import conformations_microservice as microservice
from conformations_microservice import RunPdbSimulationsRequest, OpenMmForceFields, OpenMmWaterForceFields, \
    GromacsWaterForceFields, \
    GromacsForceFields

import nolabs.api_models.conformations as api
from nolabs.api_models.experiment import TimelineResponse
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.modules.conformations.services.file_management import FileManagement
from nolabs.infrastructure.settings import Settings
from nolabs.utils import generate_uuid, utcnow


class RunSimulationsFeature:
    def __init__(self,
                 file_management: FileManagement,
                 settings: Settings):
        self._file_management = file_management
        self._settings = settings

    async def handle(self,
                     request: api.RunSimulationsRequest) -> api.RunSimulationsResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id) if request.experiment_id else ExperimentId(generate_uuid())
        experiment_name = ExperimentName(request.experiment_name)

        self._file_management.ensure_experiment_folder_exists(experiment_id)
        self._file_management.cleanup_experiment(experiment_id)
        await self._file_management.set_metadata(experiment_id, experiment_name)
        await self._file_management.set_properties(experiment_id, request)

        configuration = microservice.Configuration(
            host=self._settings.conformations_host,
        )
        with microservice.ApiClient(configuration=configuration) as _:
            client = microservice.DefaultApi(api_client=_)
            file_content = (await request.pdb_file.read()).decode('utf-8')

            timeline: List[TimelineResponse] = []

            pdb_contents = [file_content]
            fixed_file_content = self._fix_pdb(client, file_content, request, timelines=timeline)
            if fixed_file_content:
                pdb_contents.append(fixed_file_content)

            for pdb_content in pdb_contents:
                for water_ff in [GromacsWaterForceFields(wff) for wff in microservice.GromacsWaterForceFields]:
                    for ff in [GromacsForceFields(ff) for ff in microservice.GromacsForceFields]:
                        generate_gromacs_files_response = client.gen_gro_top_endpoint_gen_gro_top_post(
                            gen_gro_top_request=microservice.GenGroTopRequest(
                                ignore_missing_atoms=request.ignore_missing_atoms,
                                force_field=ff,
                                water_force_field=water_ff,
                                pdb_content=pdb_content
                            ))

                        for error in generate_gromacs_files_response.errors:
                            timeline.append(
                                TimelineResponse('Generate gromacs files error', error, created_at=utcnow())
                            )

                        if generate_gromacs_files_response.gro and generate_gromacs_files_response.top:
                            timeline.append(TimelineResponse(
                                message=f'Gromacs simulations started, force field: {ff}, water force field: {water_ff}',
                                error=None,
                                created_at=utcnow()
                            ))
                            gromacs_simulations_result = self._run_gromacs_simulations(client,
                                                                                       generate_gromacs_files_response.gro,
                                                                                       generate_gromacs_files_response.top,
                                                                                       request)

                            for error in gromacs_simulations_result.errors:
                                timeline.append(
                                    TimelineResponse(
                                        message='Gromacs simulation error',
                                        error=error,
                                        created_at=utcnow()
                                    )
                                )

                            if gromacs_simulations_result.pdb_content:
                                await self._file_management.set_properties(experiment_id=experiment_id, request=request)
                                await self._file_management.set_result(experiment_id,
                                                                       gromacs_simulations_result.pdb_content,
                                                                       timeline=timeline,
                                                                       request=request)
                                return api.RunSimulationsResponse(experiment_id=experiment_id.value,
                                                                  experiment_name=experiment_name.value,
                                                                  pdb_content=gromacs_simulations_result.pdb_content,
                                                                  timeline=timeline)

                for open_mm_water_ff in [OpenMmWaterForceFields(wff) for wff in microservice.OpenMmWaterForceFields]:
                    for open_mm_ff in [OpenMmForceFields(ff) for ff in microservice.OpenMmForceFields]:
                        timeline.append(
                            TimelineResponse(
                                message=f'Pdb simulations started, force field: {open_mm_ff}, water force field: {open_mm_water_ff}',
                                error=None,
                                created_at=utcnow()
                            )
                        )
                        pdb_simulations_result = self._run_pdb_simulations(client, pdb_content, open_mm_ff,
                                                                           open_mm_water_ff, request)

                        timeline.append(
                            TimelineResponse(
                                message=f'Pdb simulations finish',
                                error=None,
                                created_at=utcnow()
                            )
                        )

                        for error in pdb_simulations_result.errors:
                            timeline.append(
                                TimelineResponse(
                                    message='Pdb simulations error',
                                    error=error,
                                    created_at=utcnow()
                                )
                            )

                        if pdb_simulations_result.pdb_content:
                            await self._file_management.set_result(experiment_id, pdb_simulations_result.pdb_content,
                                                                   timeline,
                                                                   request)
                            return api.RunSimulationsResponse(experiment_id=experiment_id.value,
                                                              pdb_content=pdb_simulations_result.pdb_content,
                                                              timeline=timeline,
                                                              experiment_name=request.experiment_name)

            await self._file_management.set_result(experiment_id=experiment_id, pdb_content=None,
                                                   timeline=timeline, request=request)
            return api.RunSimulationsResponse(experiment_id=experiment_id.value, timeline=timeline, pdb_content=None,
                                              experiment_name=request.experiment_name)

    def _fix_pdb(self, client: microservice.DefaultApi, pdb_content: str, request: api.RunSimulationsRequest,
                 timelines: List[TimelineResponse]) -> str | None:
        timelines.append(
            TimelineResponse(
                message='Pdb fixer started',
                error=None,
                created_at=utcnow()
            )
        )
        pdb_fixer_request = microservice.RunPdbFixerRequest(
            replace_nonstandard_residues=request.replace_non_standard_residues,
            add_missing_atoms=request.add_missing_atoms,
            add_missing_hydrogens=request.add_missing_hydrogens,
            pdb_content=pdb_content
        )
        response = client.run_pdb_fixer_endpoint_run_pdb_fixer_post(run_pdb_fixer_request=pdb_fixer_request)
        timelines.append(
            TimelineResponse(
                message='Pdb fixer finished',
                error=None,
                created_at=utcnow()
            )
        )

        if not response.pdb_content:
            timelines.append(
                TimelineResponse(
                    message='Pdb fixer error',
                    error='Unexpected error in pdb fixer',
                    created_at=utcnow()
                )
            )

        for error in response.errors:
            timelines.append(
                TimelineResponse(
                    message='Pdb fixer error',
                    error=error,
                    created_at=utcnow()
                )
            )

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
