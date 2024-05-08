__all__ = [
    'router',
]

from typing import Annotated, List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, Form, UploadFile, File

from nolabs.refined.application.use_cases.ligands.api_models import *
from nolabs.refined.application.use_cases.ligands.di import LigandsControllerDependencies
from nolabs.refined.application.use_cases.ligands.use_cases import *

router = APIRouter(
    prefix='/api/v1/objects/ligands',
    tags=['Ligands'],

)


@router.post('/search',
             summary='Search ligands')
async def search_ligands(
        feature: Annotated[SearchLigandsFeature, Depends(LigandsControllerDependencies.search_ligands)],
        query: LigandSearchQuery
) -> List[LigandResponse]:
    return await feature.handle(query=query)


@router.get('/{ligand_id}',
            summary='Get ligand by id')
async def get_ligand(ligand_id: UUID, feature: Annotated[
    GetLigandFeature, Depends(LigandsControllerDependencies.get_ligand)]) -> Optional[LigandResponse]:
    return await feature.handle(ligand_id=ligand_id)


@router.post('',
             summary='Upload ligand')
async def upload_ligand(
        feature: Annotated[
            UploadLigandFeature, Depends(LigandsControllerDependencies.upload_ligand)],
        experiment_id: UUID = Form(),
        name: Optional[str] = Form(None),
        smiles: UploadFile = File(None),
        sdf: UploadFile = File(None), ) -> LigandResponse:
    return await feature.handle(request=UploadLigandRequest(
        experiment_id=experiment_id,
        name=name,
        smiles=smiles,
        sdf=sdf
    ))


@router.delete('/{ligand_id}', summary='Delete ligand')
async def delete_ligand(ligand_id: UUID, feature: Annotated[
    DeleteLigandFeature, Depends(LigandsControllerDependencies.delete_ligand)]):
    return await feature.handle(ligand_id=ligand_id)
