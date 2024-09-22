__all__ = [
    "router",
]

from typing import Annotated, List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, File, Form, UploadFile

from nolabs.application.proteins.api_models import (ProteinContentResponse,
                                                    ProteinMetadataResponse,
                                                    ProteinSearchMetadataQuery,
                                                    ProteinSearchQuery,
                                                    UpdateProteinRequest,
                                                    UploadProteinRequest,
                                                    UploadProteinResponse)
from nolabs.application.proteins.di import ProteinsControllerDependencies
from nolabs.application.proteins.use_cases import (
    DeleteProteinFeature, GetProteinFeature, GetProteinMetadataFeature,
    SearchProteinsContentFeature, SearchProteinsMetadataFeature,
    UpdateProteinFeature, UploadProteinFeature)

router = APIRouter(prefix="/api/v1/proteins", tags=["Proteins"])


@router.post("/search/content", summary="Search proteins content")
async def search_proteins(
    feature: Annotated[
        SearchProteinsContentFeature,
        Depends(ProteinsControllerDependencies.search_proteins_content),
    ],
    query: ProteinSearchQuery,
) -> List[ProteinContentResponse]:
    return await feature.handle(query=query)


@router.post("/search/metadata", summary="Search proteins metadata")
async def search_proteins(
    feature: Annotated[
        SearchProteinsMetadataFeature,
        Depends(ProteinsControllerDependencies.search_proteins_metadata),
    ],
    query: ProteinSearchMetadataQuery,
) -> List[ProteinMetadataResponse]:
    return await feature.handle(query=query)


@router.get("/{protein_id}/content", summary="Get protein content by id")
async def get_protein_content(
    protein_id: UUID,
    feature: Annotated[
        GetProteinFeature, Depends(ProteinsControllerDependencies.get_protein)
    ],
) -> Optional[ProteinContentResponse]:
    return await feature.handle(protein_id=protein_id)


@router.get("/{protein_id}/metadata", summary="Get protein metadata by id")
async def get_protein_metadata(
    protein_id: UUID,
    feature: Annotated[
        GetProteinMetadataFeature,
        Depends(ProteinsControllerDependencies.get_protein_metadata),
    ],
) -> Optional[ProteinMetadataResponse]:
    return await feature.handle(protein_id=protein_id)


@router.post("", summary="Upload protein")
async def upload_protein(
    feature: Annotated[
        UploadProteinFeature, Depends(ProteinsControllerDependencies.upload_protein)
    ],
    experiment_id: UUID = Form(),
    name: Optional[str] = Form(None),
    fasta: UploadFile = File(None),
    pdb: UploadFile = File(None),
) -> UploadProteinResponse:
    return await feature.handle(
        request=UploadProteinRequest(
            experiment_id=experiment_id, name=name, fasta=fasta, pdb=pdb
        )
    )


@router.patch("", summary="Update protein")
async def update_protein(
    feature: Annotated[
        UpdateProteinFeature, Depends(ProteinsControllerDependencies.update_protein)
    ],
    protein_id: UUID = Form(),
    name: Optional[str] = Form(None),
    fasta: UploadFile = File(None),
    pdb: UploadFile = File(None),
) -> ProteinContentResponse:
    return await feature.handle(
        request=UpdateProteinRequest(
            protein_id=protein_id, name=name, fasta=fasta, pdb=pdb
        )
    )


@router.delete("/{protein_id}", summary="Delete protein")
async def delete_protein(
    protein_id: UUID,
    feature: Annotated[
        DeleteProteinFeature, Depends(ProteinsControllerDependencies.delete_protein)
    ],
):
    return await feature.handle(protein_id=protein_id)
