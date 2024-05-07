__all__ = [
    'GetJobFeature',
    'SetupJobFeature',
    'RunJobFeature'
]

from typing import List
from uuid import UUID

import protein_design_microservice

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.use_cases.protein_design.api_models import JobResponse, SetupJobRequest
from nolabs.refined.domain.models.common import JobId, JobName, Experiment, Protein, ProteinName
from nolabs.refined.domain.models.protein_design import ProteinDesignJob
from nolabs.utils import generate_uuid


def map_job_to_response(job: ProteinDesignJob) -> JobResponse:
    return JobResponse(
        job_id=job.iid.value,
        job_name=job.name.value,
        experiment_id=job.experiment.iid.value,
        protein_id=job.protein.iid.value,
        binder_ids=[p.iid.value for p in job.binders],
        number_of_designs=job.number_of_designs,
        timesteps=job.timesteps,
        hotspots=job.hotspots
    )


class GetJobFeature:
    """
    Use case - get job information.
    """
    async def handle(self, job_id: UUID) -> JobResponse:
        try:
            job_id = JobId(job_id)
            job: ProteinDesignJob = ProteinDesignJob.objects.with_id(job_id.value)

            if not job:
                raise NoLabsException(ErrorCodes.job_not_found)

            return map_job_to_response(job)
        except Exception as e:
            if isinstance(e, NoLabsException):
                raise e

            raise NoLabsException(ErrorCodes.unknown_exception) from e


class SetupJobFeature:
    """
    Use case - create new or update existing job.
    """
    async def handle(self, request: SetupJobRequest) -> JobResponse:
        try:
            assert request

            job_id = JobId(request.job_id if request.job_id else generate_uuid())
            job_name = JobName(request.job_name if request.job_name else 'New protein design job')

            experiment = Experiment.objects.with_id(request.experiment_id)

            if not experiment:
                raise NoLabsException(ErrorCodes.experiment_not_found)

            job: ProteinDesignJob = ProteinDesignJob.objects.with_id(job_id.value)

            if not job:
                job = ProteinDesignJob(
                    id=job_id,
                    name=job_name,
                    experiment=experiment
                )

            protein = job.protein

            if not protein:
                raise NoLabsException(ErrorCodes.protein_not_found)

            job.set_input(
                protein=protein,
                contig=job.contig,
                number_of_designs=job.number_of_designs,
                hotspots=job.hotspots,
                timesteps=job.timesteps
            )
            job.save(cascade=True)

            return map_job_to_response(job)
        except Exception as e:
            print(e)
            if not isinstance(e, NoLabsException):
                raise NoLabsException(ErrorCodes.unknown_exception) from e
            raise e


class RunJobFeature:
    """
    Use case - start job.
    """
    _api: protein_design_microservice.DefaultApi

    def __init__(self, api: protein_design_microservice.DefaultApi):
        self._api = api

    async def handle(self, job_id: UUID) -> JobResponse:
        try:
            assert job_id

            job_id = JobId(job_id)
            job: ProteinDesignJob = ProteinDesignJob.objects.with_id(job_id.value)

            if not job:
                raise NoLabsException(ErrorCodes.job_not_found)

            response = self._api.run_rfdiffusion_endpoint_run_rfdiffusion_post(
                run_rfdiffusion_request=protein_design_microservice.RunRfdiffusionRequest(
                    pdb_content=job.protein.get_pdb(),
                    hotspots=job.hotspots,
                    contig=job.contig,
                    timesteps=job.timesteps,
                    number_of_designs=job.number_of_designs
                )
            )

            if response.errors and not response.pdbs_content:
                raise NoLabsException(ErrorCodes.protein_design_run_error, response.errors)

            binders: List[Protein] = []
            for i, pdb in enumerate(response.pdbs_content):
                binder = Protein.create(
                    experiment=job.experiment,
                    name=ProteinName(f'{job.protein.name.value}-binder-{str(i)}')
                )
                binder.set_pdb(pdb)
                binders.append(binder)

                job.protein.add_protein_binder(binder)
                job.save(cascade=True)
                binder.save(cascade=True)

            job.set_result(protein=job.protein, binders=binders)
            job.save(cascade=True)

            return map_job_to_response(job)
        except Exception as e:
            print(e)
            if not isinstance(e, NoLabsException):
                raise NoLabsException(ErrorCodes.unknown_exception) from e
            raise e