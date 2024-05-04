__all__ = [
    'GetProteinDesignJobFeature',
    'RunProteinDesignFeature'
]


from uuid import UUID

import protein_design_microservice

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.protein_design.api_models import GetJobResponse, JobPropertiesResponse, \
    RunProteinDesignRequest, RunProteinDesignResponse
from nolabs.refined.domain.models import JobId, JobName, Experiment, Protein, ProteinName
from nolabs.refined.domain.models.protein_design import ProteinDesignJob


class GetProteinDesignJobFeature:
    async def handle(self, experiment_id: UUID, job_id: UUID) -> GetJobResponse:
        assert job_id
        assert experiment_id

        experiment = Experiment.objects.with_id(id=experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        job_id = JobId(value=job_id)

        job: ProteinDesignJob = ProteinDesignJob.objects.with_id(id=job_id.value)

        if not job:
            job = ProteinDesignJob(
                id=job_id,
                name=JobName('New protein design job'),
                experiment=experiment
            )
            job.save()
            return GetJobResponse(
                job_id=job.iid.value,
                job_name=job.name.value,
                pdb_files=[],
                properties=JobPropertiesResponse.default()
            )

        return GetJobResponse(
            job_id=job.iid.value,
            job_name=job.name.value,
            pdb_files=[binder.get_pdb() for binder in job.binders],
            properties=JobPropertiesResponse(
                pdb_file=job.protein.get_pdb(),
                contig=job.contig,
                number_of_designs=job.number_of_designs,
                hotspots=job.hotspots,
                timesteps=job.timesteps,
                pdb_file_name=job.protein.name.pdb_name
            )
        )


class RunProteinDesignFeature:
    def __init__(self, api: protein_design_microservice.DefaultApi):
        self._api = api

    async def handle(self,
                     request: RunProteinDesignRequest) -> RunProteinDesignResponse:
        assert request

        job_id = JobId(request.job_id)
        job: ProteinDesignJob = ProteinDesignJob.objects.with_id(id=job_id.value)

        protein = Protein.create(
            experiment=job.experiment,
            name=ProteinName(value=request.pdb_file.filename),
            pdb_content=await request.pdb_file.read()
        )

        job.set_input(
            protein=protein,
            contig=request.contig,
            number_of_designs=request.number_of_designs,
            hotspots=request.hotspots,
            timesteps=request.timesteps
        )

        job.save(cascade=True)

        response = self._api.run_rfdiffusion_endpoint_run_rfdiffusion_post(
            run_rfdiffusion_request=protein_design_microservice.RunRfdiffusionRequest(
                pdb_content=job.protein.get_pdb(),
                hotspots=request.hotspots,
                contig=request.contig,
                timesteps=request.timesteps,
                number_of_designs=request.number_of_designs
            )
        )

        if response.errors and not response.pdbs_content:
            raise NoLabsException(ErrorCodes.protein_design_run_error, response.errors)

        for i, pdb in enumerate(response.pdbs_content):
            binder = Protein.create(
                experiment=job.experiment,
                name=ProteinName(f'{protein.name.value}-binder-{str(i)}'),
                pdb_content=pdb
            )
            job.set_result(protein, binder)
            protein.set_protein_binder(binder)
            binder.set_protein_binder(protein)

        job.save(cascade=True)

        return RunProteinDesignResponse(
            job_id=job.iid.value,
            pdb_files=response.pdbs_content
        )
