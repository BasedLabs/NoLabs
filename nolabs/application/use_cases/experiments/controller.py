__all__ = ["router"]

from typing import Annotated, List
from uuid import UUID

from fastapi import APIRouter, Depends

from nolabs.application.use_cases.experiments.api_models import (
    ExperimentMetadataResponse, UpdateExperimentRequest)
from nolabs.application.use_cases.experiments.di import ExperimentsDependencies
from nolabs.application.use_cases.experiments.use_cases import (
    CreateExperimentFeature, DeleteExperimentFeature,
    GetExperimentsMetadataFeature, UpdateExperimentFeature)

router = APIRouter(prefix="/api/v1/experiments", tags=["Experiments"])


@router.get("/all", summary="Get all experiments")
async def experiments(
    feature: Annotated[
        GetExperimentsMetadataFeature, Depends(ExperimentsDependencies.experiments)
    ]
) -> List[ExperimentMetadataResponse]:
    return feature.handle()


@router.delete("/{experiment_id}", summary="Delete experiment")
async def delete_experiment(
    experiment_id: UUID,
    feature: Annotated[
        DeleteExperimentFeature, Depends(ExperimentsDependencies.delete_experiment)
    ],
):
    return await feature.handle(experiment_id)


@router.patch("/", summary="Update experiment")
async def update_experiment(
    request: UpdateExperimentRequest,
    feature: Annotated[
        UpdateExperimentFeature, Depends(ExperimentsDependencies.update_experiment)
    ],
):
    return feature.handle(request)


@router.post("/", summary="Create experiment")
async def create_experiment(
    feature: Annotated[
        CreateExperimentFeature, Depends(ExperimentsDependencies.create_experiment)
    ]
) -> ExperimentMetadataResponse:
    return feature.handle()
