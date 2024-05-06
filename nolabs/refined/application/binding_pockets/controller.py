__all__ = [
    'router',
]

from typing import List, Annotated, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, Form, File, UploadFile

from nolabs.refined.application.binding_pockets.api_models import PredictBindingPocketsRequest, BindingPocketsResponse
from nolabs.refined.application.binding_pockets.use_cases import PredictBindingPocketsFeature, GetBindingPocketsFeature, \
    SetManualBindingPocketsFeature
from nolabs.refined.application.binding_pockets.di import BindingPocketsDependencies

router = APIRouter(
    prefix='/api/v1/binding-pockets',
    tags=['Binding pockets'],

)


@router.post('/jobs/start',
             summary='Start binding pockets prediction job')
async def predict(
        feature: Annotated[
            PredictBindingPocketsFeature, Depends(BindingPocketsDependencies.predict_binding_pockets)],
        request: PredictBindingPocketsRequest
) -> BindingPocketsResponse:
    return await feature.handle(request)


@router.get('/', summary='Get binding pockets')
async def get_binding_pockets(protein_id: UUID, feature: Annotated[
    GetBindingPocketsFeature, Depends(BindingPocketsDependencies.get_binding_pockets)]) -> BindingPocketsResponse:
    return await feature.handle(protein_id=protein_id)


@router.post('/', summary='Set binding pockets')
async def get_binding_pockets(protein_id: UUID, binding_pockets: List[int], feature: Annotated[
    SetManualBindingPocketsFeature, Depends(BindingPocketsDependencies.set_binding_pockets)]) -> BindingPocketsResponse:
    return await feature.handle(protein_id=protein_id, binding_pockets=binding_pockets)
