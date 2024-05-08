__all__ = [
    'router',
]

from typing import Annotated, List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, Form, UploadFile, File

from nolabs.refined.application.use_cases.protein_design.api_models import SetupJobRequest, JobResponse
from nolabs.refined.application.use_cases.protein_design.di import ProteinDesignDependencies
from nolabs.refined.application.use_cases.protein_design.use_cases import RunJobFeature, GetJobFeature, SetupJobFeature
from nolabs.refined.application.use_cases.proteins.api_models import ProteinSearchQuery, ProteinResponse, \
    UploadProteinRequest
from nolabs.refined.application.use_cases.proteins.di import ProteinsControllerDependencies
from nolabs.refined.application.use_cases.proteins.use_cases import SearchProteinsFeature, GetProteinFeature, \
    UploadProteinFeature, DeleteProteinFeature

router = APIRouter(
    prefix='/api/v1/objects/proteins',
    tags=['Proteins'],

)


@router.post('/search',
             summary='Search proteins')
async def search_proteins(
        feature: Annotated[SearchProteinsFeature, Depends(ProteinsControllerDependencies.search_proteins)],
        query: ProteinSearchQuery
) -> List[ProteinResponse]:
    return await feature.handle(query=query)


@router.get('/{protein_id}',
            summary='Get protein by id')
async def get_protein(protein_id: UUID, feature: Annotated[
    GetProteinFeature, Depends(ProteinsControllerDependencies.get_protein)]) -> Optional[ProteinResponse]:
    return await feature.handle(protein_id=protein_id)


@router.post('',
             summary='Upload protein')
async def upload_protein(
        feature: Annotated[
            UploadProteinFeature, Depends(ProteinsControllerDependencies.upload_protein)],
        experiment_id: UUID = Form(),
        name: Optional[str] = Form(None),
        fasta: UploadFile = File(None),
        pdb: UploadFile = File(None), ) -> ProteinResponse:
    return await feature.handle(request=UploadProteinRequest(
        experiment_id=experiment_id,
        name=name,
        fasta=fasta,
        pdb=pdb
    ))


@router.delete('/{protein_id}', summary='Delete protein')
async def delete_protein(protein_id: UUID, feature: Annotated[
    DeleteProteinFeature, Depends(ProteinsControllerDependencies.delete_protein)]):
    return await feature.handle(protein_id=protein_id)
