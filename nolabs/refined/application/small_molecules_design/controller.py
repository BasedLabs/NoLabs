from typing import List, Annotated, Optional
from uuid import UUID

from fastapi import WebSocket, APIRouter, UploadFile, File, Depends

from nolabs.refined.application.small_molecules_design.api_models import GetJobStatusResponse, LogsResponse, \
    SmilesResponse, JobResponse, SetupJobRequest
from nolabs.refined.application.small_molecules_design.di import SmallMoleculesDesignDependencies
from nolabs.refined.application.small_molecules_design.use_cases import GetJobStatus, GetJobFeature, StopJobFeature, \
    DeleteJobFeature, GetJobLogsFeature, \
    GetJobSmilesFeature, RunLearningStageJobFeature, RunSamplingStageJobFeature, SetupJobFeature

logs_websocket: WebSocket | None

router = APIRouter(
    prefix='/api/v1/small-molecules-design',
    tags=['small-molecules-design']
)


@router.get('/jobs/{job_id}/status')
async def get_job_status(
        feature: Annotated[GetJobStatus, Depends(SmallMoleculesDesignDependencies.get_job_status)],
        job_id: UUID
) -> GetJobStatusResponse:
    return await feature.handle(job_id=job_id)


@router.post("/jobs/{job_id}/run/learning")
async def run_learning_stage_job(
        feature: Annotated[RunLearningStageJobFeature, Depends(SmallMoleculesDesignDependencies.run_learning)],
        job_id: UUID):
    return await feature.handle(job_id=job_id)


@router.post("/jobs/{job_id}/run/sampling")
async def run_sampling_stage_job(
        feature: Annotated[RunSamplingStageJobFeature, Depends(SmallMoleculesDesignDependencies.run_sampling)],
        job_id: UUID):
    return await feature.handle(job_id=job_id)


@router.get('/jobs/{job_id}')
async def get_job(feature: Annotated[GetJobFeature, Depends(SmallMoleculesDesignDependencies.get_job)],
                         job_id: UUID) -> Optional[JobResponse]:
    return await feature.handle(job_id)


@router.post('/jobs')
async def setup_job(
        feature: Annotated[SetupJobFeature, Depends(SmallMoleculesDesignDependencies.setup_job)],
        request: SetupJobRequest) -> JobResponse:
    return await feature.handle(request=request)


@router.post("/jobs/{job_id}/stop")
async def stop_job(
        feature: Annotated[StopJobFeature, Depends(SmallMoleculesDesignDependencies.stop_job)],
        job_id: UUID):
    return await feature.handle(job_id=job_id)


@router.delete("/jobs/{job_id}")
async def delete_job(
        feature: Annotated[DeleteJobFeature, Depends(SmallMoleculesDesignDependencies.delete_job)],
        job_id: UUID):
    return await feature.handle(job_id=job_id)


@router.get('/jobs/{job_id}/logs')
async def get_jog_logs(feature: Annotated[GetJobLogsFeature, Depends(SmallMoleculesDesignDependencies.get_job_logs)],
               job_id: UUID) -> LogsResponse:
    return await feature.handle(job_id=job_id)


@router.get('/jobs/{job_id}/smiles')
async def get_job_smiles(feature: Annotated[GetJobSmilesFeature, Depends(SmallMoleculesDesignDependencies.get_job_smiles)],
                 job_id: UUID) -> List[SmilesResponse]:
    return await feature.handle(job_id=job_id)

