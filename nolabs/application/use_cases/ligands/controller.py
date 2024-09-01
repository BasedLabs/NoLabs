__all__ = [
    'router',
]

from typing import Annotated, List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, Form, UploadFile, File

from nolabs.application.use_cases.ligands.api_models import *
from nolabs.application.use_cases.ligands.api_models import UpdateLigandRequest, UploadLigandResponse
from nolabs.application.use_cases.ligands.di import LigandsControllerDependencies
from nolabs.application.use_cases.ligands.use_cases import *
from nolabs.application.use_cases.ligands.use_cases import UpdateLigandFeature

router = APIRouter(
    prefix='/api/v1/objects/ligands',
    tags=['Ligands'],

)


@router.post('/search/content',
             summary='Search ligands content')
async def search_ligands_content(
        feature: Annotated[SearchLigandsContentFeature, Depends(LigandsControllerDependencies.search_ligands_content)],
        query: LigandSearchContentQuery
) -> List[LigandContentResponse]:
    return await feature.handle(query=query)


@router.post('/search/metadata',
             summary='Search ligands metadata')
async def search_ligands_metadata(
        feature: Annotated[SearchLigandsMetadataFeature, Depends(LigandsControllerDependencies.search_ligands_metadata)],
        query: LigandSearchMetadataQuery
) -> List[LigandMetadataResponse]:
    return await feature.handle(query=query)


@router.get('/{ligand_id}/content',
            summary='Get ligand content by id')
async def get_ligand_content(ligand_id: UUID, feature: Annotated[
    GetLigandFeature, Depends(LigandsControllerDependencies.get_ligand)]) -> Optional[LigandContentResponse]:
    return await feature.handle(ligand_id=ligand_id)


@router.get('/{ligand_id}/content',
            summary='Get ligand content by id')
async def get_ligand_content(ligand_id: UUID, feature: Annotated[
    GetLigandFeature, Depends(LigandsControllerDependencies.get_ligand)]) -> Optional[LigandContentResponse]:
    return await feature.handle(ligand_id=ligand_id)


@router.post('',
             summary='Upload ligand')
async def upload_ligand(
        feature: Annotated[
            UploadLigandFeature, Depends(LigandsControllerDependencies.upload_ligand)],
        experiment_id: UUID = Form(),
        name: Optional[str] = Form(None),
        smiles: UploadFile = File(None),
        sdf: UploadFile = File(None), ) -> UploadLigandResponse:
    return await feature.handle(request=UploadLigandRequest(
        experiment_id=experiment_id,
        name=name,
        smiles=smiles,
        sdf=sdf
    ))


@router.patch('',
             summary='Update ligand')
async def update_ligand(
        feature: Annotated[
            UpdateLigandFeature, Depends(LigandsControllerDependencies.update_ligand)],
        ligand_id: UUID = Form(),
        name: Optional[str] = Form(None),
        smiles: UploadFile = File(None),
        sdf: UploadFile = File(None), ) -> LigandContentResponse:
    return await feature.handle(request=UpdateLigandRequest(
        ligand_id=ligand_id,
        name=name,
        smiles=smiles,
        sdf=sdf
    ))


@router.delete('/{ligand_id}', summary='Delete ligand')
async def delete_ligand(ligand_id: UUID, feature: Annotated[
    DeleteLigandFeature, Depends(LigandsControllerDependencies.delete_ligand)]):
    return await feature.handle(ligand_id=ligand_id)
