import datetime
from typing import List
from uuid import UUID

import conformations_microservice

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.conformations.api_models import GetJobResponse, JobPropertiesResponse, \
    RunSimulationsRequest, RunSimulationsResponse, IntegratorsRequest
from nolabs.refined.application.experiments.api_models import TimelineResponse
from nolabs.refined.domain.models import JobId, JobName, Experiment, Protein, ProteinName
from nolabs.refined.domain.models.conformations import ConformationsJob, Integrator, ConformationsTimeline
from nolabs.utils import generate_uuid


def map_timeline(timeline: TimelineResponse) -> ConformationsTimeline:
    return ConformationsTimeline(
        message=timeline.message,
        error=timeline.error,
        created_at=timeline.created_at
    )


class GetConformationsJobFeature:
    async def handle(self, job_id: UUID) -> GetJobResponse:
        assert id

        job_id = JobId(value=job_id)
        job: ConformationsJob = ConformationsJob.objects.with_id(id=job_id.value)

        if job:
            raise NoLabsException(ErrorCodes.job_not_found)

        return GetJobResponse(
            job_id=job.iid.value,
            job_name=job.name.value,
            pdb_file=job.protein.get_md(),
            timeline=[TimelineResponse(
                message=timeline.message,
                error=timeline.error,
                created_at=timeline.created_at
            ) for timeline in job.timeline],
            properties=JobPropertiesResponse(
                pdb_file=job.protein.get_pdb(),
                pdb_file_name=job.protein.name.value,
                total_frames=job.total_frames,
                temperature_k=job.temperature_k,
                take_frame_every=job.take_frame_every,
                step_size=job.step_size,
                replace_non_standard_residues=job.replace_non_standard_residues,
                add_missing_atoms=job.add_missing_atoms,
                add_missing_hydrogens=job.add_missing_hydrogens,
                friction_coeff=job.friction_coeff,
                ignore_missing_atoms=job.ignore_missing_atoms,
                integrator=IntegratorsRequest(job.integrator)
            )
        )


class RunConformationsJobFeature:
    def __init__(self,
                 api: conformations_microservice.DefaultApi):
        self._api = api

    async def handle(self,
                     request: RunSimulationsRequest) -> RunSimulationsResponse:
        assert request

        if not request.job_id:
            request.job_id = generate_uuid()

        experiment = Experiment.objects.with_id(id=request.experiment_id)
        job = ConformationsJob.objects.with_id(id=request.job_id)

        if not job:
            job = ConformationsJob(
                id=JobId(request.job_id),
                name=JobName('New molecular dynamics job'),
                experiment=experiment
            )

        protein_content = await request.pdb_file.read()
        protein_name = ProteinName(request.pdb_file.filename)
        protein = Protein.create(
            experiment=experiment,
            name=protein_name,
            pdb_content=protein_content
        )
        job.set_input(
            protein=protein,
            integrator=Integrator(request.integrator.value),
            total_frames=request.total_frames,
            temperature_k=request.temperature_k,
            take_frame_every=request.take_frame_every,
            step_size=request.step_size,
            replace_non_standard_residues=request.replace_non_standard_residues,
            add_missing_atoms=request.add_missing_atoms,
            add_missing_hydrogens=request.add_missing_hydrogens,
            friction_coeff=request.friction_coeff,
            ignore_missing_atoms=request.ignore_missing_atoms
        )

        timeline: List[TimelineResponse] = []

        pdb_contents = [protein_content]
        fixed_file_content = self._fix_pdb(protein.get_pdb(), request, timelines=timeline)
        if fixed_file_content:
            pdb_contents.append(fixed_file_content)

        job.clear_result()

        for pdb_content in pdb_contents:
            for water_ff in [conformations_microservice.GromacsWaterForceFields(wff) for wff in
                             conformations_microservice.GromacsWaterForceFields]:
                for ff in [conformations_microservice.GromacsForceFields(ff) for ff in
                           conformations_microservice.GromacsForceFields]:
                    generate_gromacs_files_response = self._api.gen_gro_top_endpoint_gen_gro_top_post(
                        gen_gro_top_request=conformations_microservice.GenGroTopRequest(
                            ignore_missing_atoms=request.ignore_missing_atoms,
                            force_field=ff,
                            water_force_field=water_ff,
                            pdb_content=pdb_content
                        ))

                    for error in generate_gromacs_files_response.errors:
                        tl_entry = TimelineResponse('Generate gromacs files error', error,
                                                    created_at=datetime.datetime.utcnow())
                        timeline.append(tl_entry)
                        job.append_timeline(protein=protein, timeline=map_timeline(tl_entry))
                        job.save()

                    if generate_gromacs_files_response.gro and generate_gromacs_files_response.top:
                        tl_entry = TimelineResponse(
                            message=f'Gromacs simulations started, force field: {ff}, water force field: {water_ff}',
                            error=None,
                            created_at=datetime.datetime.utcnow()
                        )
                        timeline.append(tl_entry)
                        job.append_timeline(protein=protein, timeline=map_timeline(tl_entry))
                        job.save()
                        gromacs_simulations_result = self._run_gromacs_simulations(generate_gromacs_files_response.gro,
                                                                                   generate_gromacs_files_response.top,
                                                                                   request)

                        for error in gromacs_simulations_result.errors:
                            tl_entry = TimelineResponse(
                                message='Gromacs simulation error',
                                error=error,
                                created_at=datetime.datetime.utcnow()
                            )
                            timeline.append(tl_entry)
                            job.append_timeline(protein=protein, timeline=map_timeline(tl_entry))
                            job.save()

                        if gromacs_simulations_result.pdb_content:
                            job.set_result(protein=protein,
                                           md_content=gromacs_simulations_result.pdb_content)
                            protein.set_protein(gromacs_simulations_result.pdb_content)
                            job.save(cascade=True)
                            protein.save(cascade=True)
                            return RunSimulationsResponse(job_id=job.iid.value,
                                                          pdb_content=gromacs_simulations_result.pdb_content,
                                                          timeline=timeline)

            for open_mm_water_ff in [conformations_microservice.OpenMmWaterForceFields(wff) for wff in
                                     conformations_microservice.OpenMmWaterForceFields]:
                for open_mm_ff in [conformations_microservice.OpenMmForceFields(ff) for ff in
                                   conformations_microservice.OpenMmForceFields]:
                    tl_entry = TimelineResponse(
                        message=f'Pdb simulations started, force field: {open_mm_ff}, water force field: {open_mm_water_ff}',
                        error=None,
                        created_at=datetime.datetime.utcnow()
                    )
                    timeline.append(tl_entry)
                    job.append_timeline(protein=protein, timeline=map_timeline(tl_entry))
                    pdb_simulations_result = self._run_pdb_simulations(pdb_content, open_mm_ff,
                                                                       open_mm_water_ff, request)

                    tl_entry = TimelineResponse(
                        message=f'Pdb simulations finish',
                        error=None,
                        created_at=datetime.datetime.utcnow()
                    )
                    timeline.append(tl_entry)
                    job.append_timeline(protein=protein, timeline=map_timeline(tl_entry))
                    job.save()

                    for error in pdb_simulations_result.errors:
                        tl_entry = TimelineResponse(
                            message='Pdb simulations error',
                            error=error,
                            created_at=datetime.datetime.utcnow()
                        )
                        timeline.append(tl_entry)
                        job.append_timeline(protein=protein, timeline=map_timeline(tl_entry))
                        job.save()

                    if pdb_simulations_result.pdb_content:
                        job.set_result(protein=protein,
                                       md_content=pdb_simulations_result.pdb_content)
                        protein.set_md_content(pdb_simulations_result.pdb_content)
                        job.save(cascade=True)
                        protein.save()
                        return RunSimulationsResponse(job_id=job.iid.value,
                                                      pdb_content=pdb_simulations_result.pdb_content,
                                                      timeline=timeline)

        return RunSimulationsResponse(job_id=job.iid.value, timeline=timeline, pdb_content=None)

    def _fix_pdb(self, pdb_content: str, request: RunSimulationsRequest,
                 timelines: List[TimelineResponse]) -> str | None:
        timelines.append(
            TimelineResponse(
                message='Pdb fixer started',
                error=None,
                created_at=datetime.datetime.utcnow()
            )
        )
        pdb_fixer_request = conformations_microservice.RunPdbFixerRequest(
            replace_nonstandard_residues=request.replace_non_standard_residues,
            add_missing_atoms=request.add_missing_atoms,
            add_missing_hydrogens=request.add_missing_hydrogens,
            pdb_content=pdb_content
        )
        response = self._api.run_pdb_fixer_endpoint_run_pdb_fixer_post(run_pdb_fixer_request=pdb_fixer_request)
        timelines.append(
            TimelineResponse(
                message='Pdb fixer finished',
                error=None,
                created_at=datetime.datetime.utcnow()
            )
        )

        if not response.pdb_content:
            timelines.append(
                TimelineResponse(
                    message='Pdb fixer error',
                    error='Unexpected error in pdb fixer',
                    created_at=datetime.datetime.utcnow()
                )
            )

        for error in response.errors:
            timelines.append(
                TimelineResponse(
                    message='Pdb fixer error',
                    error=error,
                    created_at=datetime.datetime.utcnow()
                )
            )

        return response.pdb_content

    def _run_gromacs_simulations(self, gro: str, top: str,
                                 request: RunSimulationsRequest) -> conformations_microservice.RunSimulationsResponse:
        run_gromacs_simulations_request = conformations_microservice.RunGromacsSimulationsRequest(
            temperature_k=request.temperature_k,
            friction_coeff=request.friction_coeff,
            step_size=request.step_size,
            integrator=microservice.Integrators(request.integrator.value),  # type: ignore
            take_frame_every=request.take_frame_every,
            total_frames=request.total_frames,
            top=top,
            gro=gro
        )
        return self._api.run_gromacs_simulations_endpoint_run_gromacs_simulations_post(
            run_gromacs_simulations_request=run_gromacs_simulations_request)

    def _run_pdb_simulations(self,
                             pdb_content: str,
                             ff: conformations_microservice.OpenMmForceFields,
                             water_ff: conformations_microservice.OpenMmWaterForceFields,
                             request: RunSimulationsRequest) -> conformations_microservice.RunSimulationsResponse:
        run_pdb_simulations_request = conformations_microservice.RunPdbSimulationsRequest(
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
        return self._api.run_pdb_simulations_endpoint_run_pdb_simulations_post(
            run_pdb_simulations_request=run_pdb_simulations_request)
