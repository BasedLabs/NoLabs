from typing import List, Annotated, Optional
from uuid import UUID

from fastapi import WebSocket, APIRouter, UploadFile, File, Depends

from nolabs.modules.small_molecules_design.get_smiles import GetSmilesFeature
from nolabs.refined.application.small_molecules_design.api_models import GetJobStatusResponse, SamplingSizeRequest, \
    GetJobResponse, JobPropertiesRequest, LogsResponse, SmilesResponse
from nolabs.refined.application.small_molecules_design.di import SmallMoleculesDesignDependencies
from nolabs.refined.application.small_molecules_design.use_cases import GetJobStatus, StartLearningJobFeature, \
    StartSamplingJobFeature, GetJobFeature, SavePropertiesFeature, StopJobFeature, DeleteJobFeature, GetJobLogsFeature, \
    GetJobSmilesFeature

logs_websocket: WebSocket | None

router = APIRouter(
    prefix='/api/v1/small-molecules-design',
    tags=['small-molecules-design']
)


@router.get('/jobs/{job_id}/status')
async def status(
        feature: Annotated[GetJobStatus, Depends(SmallMoleculesDesignDependencies.get_job_status)],
        job_id: UUID
) -> GetJobStatusResponse:
    return await feature.handle(job_id=job_id)


@router.post("/jobs/{job_id}/learning")
async def learning(
        feature: Annotated[StartLearningJobFeature, Depends(SmallMoleculesDesignDependencies.start_learning)],
        job_id: UUID):
    return await feature.handle(job_id=job_id)


@router.post("/jobs/{job_id}/sampling")
async def sampling(
        feature: Annotated[StartSamplingJobFeature, Depends(SmallMoleculesDesignDependencies.start_sampling)],
        job_id: UUID,
        request: SamplingSizeRequest):
    return await feature.handle(job_id=job_id, request=request)


@router.get('/jobs/{job_id}')
async def get_experiment(feature: Annotated[GetJobFeature, Depends(SmallMoleculesDesignDependencies.get_job)],
                         job_id: UUID) -> Optional[GetJobResponse]:
    return await feature.handle(job_id)


@router.post('/jobs/props')
async def save_properties(
        feature: Annotated[SavePropertiesFeature, Depends(SmallMoleculesDesignDependencies.save_properties)],
        job_id: UUID,
        experiment_id: UUID,
        center_x: float,
        center_y: float, center_z: float, size_x: float,
        size_y: float, size_z: float,
        epochs: int = 50,
        batch_size: int = 128,
        minscore: float = 0.4,
        pdb_file: UploadFile = File()):
    await feature.handle(job_id=job_id,
                         experiment_id=experiment_id,
                         request=JobPropertiesRequest(
                             pdb_file=pdb_file,
                             center_x=center_x,
                             center_y=center_y,
                             center_z=center_z,
                             size_x=size_x,
                             size_y=size_y,
                             size_z=size_z,
                             batch_size=batch_size,
                             minscore=minscore,
                             epochs=epochs
                         ))


@router.post("/jobs/{job_id}/stop")
async def stop(
        feature: Annotated[StopJobFeature, Depends(SmallMoleculesDesignDependencies.stop_job)],
        job_id: UUID):
    return await feature.handle(job_id=job_id)


@router.delete("/jobs/{job_id}")
async def delete(
        feature: Annotated[DeleteJobFeature, Depends(SmallMoleculesDesignDependencies.delete_job)],
        job_id: UUID):
    return await feature.handle(job_id=job_id)


@router.get('/jobs/{job_id}/logs')
async def logs(feature: Annotated[GetJobLogsFeature, Depends(SmallMoleculesDesignDependencies.get_job_logs)],
               job_id: UUID) -> LogsResponse:
    return await feature.handle(job_id=job_id)


@router.get('/jobs/{job_id}/smiles')
async def smiles(feature: Annotated[GetJobSmilesFeature, Depends(SmallMoleculesDesignDependencies.get_job_smiles)],
                 job_id: UUID) -> List[SmilesResponse]:
    return await feature.handle(job_id=job_id)

