__all__ = [
    'DeleteJobFeature',
    'GetJobStatus',
    'GetJobFeature',
    'GetJobLogsFeature',
    'GetJobSmilesFeature',
    'SavePropertiesFeature',
    'StartLearningJobFeature',
    'StartSamplingJobFeature',
    'StopJobFeature'
]

from typing import Optional, List
from uuid import UUID

import reinvent_microservice

from nolabs.refined.application.small_molecules_design.api_models import GetJobStatusResponse, \
    GetJobResponse, JobPropertiesResponse, LogsResponse, SmilesResponse, JobPropertiesRequest, SamplingSizeRequest
from nolabs.refined.domain.models import JobId, Job, JobName, Experiment, Protein, ProteinName
from nolabs.refined.domain.models.small_molecules_design import SmallMoleculesDesignJob


class DeleteJobFeature:
    def __init__(self, api: reinvent_microservice.ReinventApi):
        self._api = api

    async def handle(self, job_id: UUID):
        self._api.delete_api_reinvent_config_id_delete(config_id=job_id)
        job_id = JobId(job_id)
        job: Job = Job.objects.with_id(job_id.value)
        if not job:
            return
        job.delete()


class GetJobStatus:
    def __init__(self, api: reinvent_microservice.ReinventApi):
        self._api = api

    async def handle(self, job_id: UUID) -> GetJobStatusResponse:
        if not Job.objects(id=job_id):
            return GetJobStatusResponse(running=False, sampling_allowed=False)

        config_result = self._api.get_config_api_reinvent_reinvent_config_id_get(config_id=job_id)

        if not config_result:
            return GetJobStatusResponse(running=False, sampling_allowed=False)

        config = config_result.actual_instance

        return GetJobStatusResponse(running=config.running, sampling_allowed=config.sampling_allowed)


class GetJobFeature:
    def __init__(self, api: reinvent_microservice.ReinventApi):
        self._api = api

    async def handle(self, job_id: UUID) -> Optional[GetJobResponse]:
        if not Job.objects(id=job_id):
            return None

        job_id = JobId(job_id)
        job: SmallMoleculesDesignJob = SmallMoleculesDesignJob.objects.with_id(id=job_id.value)

        properties = JobPropertiesResponse(
            pdb_file=job.protein.get_pdb() if job.protein else None,
            pdb_file_name=job.protein.name.value if job.protein else None,
            center_x=job.center_x,
            center_y=job.center_y,
            center_z=job.center_z,
            size_x=job.size_x,
            size_y=job.size_y,
            size_z=job.size_z,
            batch_size=job.batch_size,
            minscore=job.minscore,
            epochs=job.epochs
        )

        config_result = self._api.get_config_api_reinvent_reinvent_config_id_get(config_id=job_id.value)

        if not config_result:
            return GetJobResponse(
                job_id=job_id.value,
                job_name=job.name.value,
                created_at=job.created_at,
                status=GetJobStatusResponse(running=False, sampling_allowed=False),
                properties=properties
            )

        config = config_result.actual_instance

        running = config.running if config_result else False
        sampling_allowed = config.sampling_allowed if config_result else False

        return GetJobResponse(
            job_id=job_id.value,
            job_name=job.name.value,
            created_at=job.created_at,
            status=GetJobStatusResponse(running=running, sampling_allowed=sampling_allowed),
            properties=properties
        )


class GetJobLogsFeature:
    def __init__(self, api: reinvent_microservice.ReinventApi):
        self._api = api

    async def handle(self, job_id: UUID) -> LogsResponse:
        if not Job.objects(id=job_id):
            return LogsResponse(
                output='',
                docking_output='',
                errors=''
            )

        response = self._api.logs_api_reinvent_config_id_logs_get(config_id=job_id)

        if not response:
            return LogsResponse(
                output='Nothing',
                docking_output='Nothing',
                errors='Nothing'
            )

        instance = response.actual_instance

        return LogsResponse(output=instance.output,
                            docking_output=instance.docking_output,
                            errors=instance.errors)


class GetJobSmilesFeature:
    def __init__(self, api: reinvent_microservice.ReinventApi):
        self._api = api

    async def handle(self, job_id: UUID) -> List[SmilesResponse]:
        job = Job.objects.with_id(id=job_id)
        if not job:
            return []

        response = self._api.smiles_api_reinvent_config_id_smiles_get(config_id=job_id)

        return [SmilesResponse(
            smiles=s.smiles,
            drug_likeness=s.drug_likeness,
            score=s.score,
            stage=s.stage,
            created_at=job.created_at
        ) for s in response.smiles]


class SavePropertiesFeature:
    def __init__(self, api: reinvent_microservice.ReinventApi):
        self._api = api

    async def handle(self, experiment_id: UUID, job_id: UUID, request: JobPropertiesRequest):
        experiment = Experiment.objects.with_id(id=experiment_id)
        job = SmallMoleculesDesignJob.objects.with_id(id=job_id)
        if not job:
            job = SmallMoleculesDesignJob(
                id=JobId(job_id),
                name=JobName('New small molecules design job'),
                experiment=experiment
            )
        pdb_file = await request.pdb_file.read()

        protein = Protein.create(
            experiment=experiment,
            name=ProteinName(request.pdb_file.filename),
            pdb_content=pdb_file
        )

        job.set_inputs(
            protein=protein,
            center_x=request.center_x,
            center_y=request.center_y,
            center_z=request.center_z,
            size_x=request.size_x,
            size_y=request.size_y,
            size_z=request.size_z,
            batch_size=request.batch_size,
            minscore=request.minscore,
            epochs=request.epochs,
        )

        job.save()

        tmp_pdb_path = 'tmp.pdb'
        open(tmp_pdb_path, 'w').write(protein.get_pdb())

        self._api.save_params_api_reinvent_config_id_params_post(
            config_id=experiment_id,
            name=job.name,
            pdb_file=tmp_pdb_path,
            center_x=request.center_x,
            center_y=request.center_y,
            center_z=request.center_z,
            size_x=request.size_x,
            size_y=request.size_y,
            size_z=request.size_z,
            batch_size=request.batch_size,
            minscore=request.minscore,
            epochs=request.epochs
        )


class StartLearningJobFeature:
    def __init__(self, api: reinvent_microservice.ReinventApi):
        self._api = api

    async def handle(self, job_id: UUID):
        self._api.learning_api_reinvent_config_id_start_learning_post(config_id=job_id)


class StartSamplingJobFeature:
    def __init__(self, api: reinvent_microservice.ReinventApi):
        self._api = api

    async def handle(self, job_id: UUID, request: SamplingSizeRequest):
        self._api.sampling_api_reinvent_config_id_start_sampling_post(config_id=job_id,
                                                                      sampling_size_request=reinvent_microservice.SamplingSizeRequest(
                                                                          number_of_molecules_to_design=request.number))


class StopJobFeature:
    def __init__(self, api: reinvent_microservice.ReinventApi):
        self._api = api

    async def handle(self, job_id: UUID):
        self._api.stop_api_reinvent_config_id_jobs_stop_post(config_id=job_id)
