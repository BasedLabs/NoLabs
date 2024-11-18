__all__ = [
    "router",
]

from typing import Annotated, List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, File, Form, UploadFile

from nolabs.application.ligands.api_models import (
    LigandContentResponse,
    LigandMetadataResponse,
    LigandSearchContentQuery,
    LigandSearchMetadataQuery,
    UpdateLigandRequest,
    UploadLigandRequest,
    UploadLigandResponse,
)
from nolabs.application.ligands.use_cases import (
    DeleteLigandFeature,
    GetLigandFeature,
    SearchLigandsContentFeature,
    SearchLigandsMetadataFeature,
    UpdateLigandFeature,
    UploadLigandFeature,
)

router = APIRouter(
    prefix="/api/v1/objects/ligands",
    tags=["Ligands"],
)


@router.post("/search/content", summary="Search ligands content")
async def search_ligands_content(query: LigandSearchContentQuery, ) -> List[LigandContentResponse]:
    return await SearchLigandsContentFeature().handle(query=query)


@router.post("/search/metadata", summary="Search ligands metadata")
async def search_ligands_metadata(query: LigandSearchMetadataQuery, ) -> List[LigandMetadataResponse]:
    return await SearchLigandsMetadataFeature().handle(query=query)


@router.get("/{ligand_id}/content", summary="Get ligand content by id")
async def get_ligand_content(ligand_id: UUID) -> Optional[LigandContentResponse]:
    return await GetLigandFeature().handle(ligand_id=ligand_id)


@router.get("/{ligand_id}/content", summary="Get ligand content by id")
async def get_ligand_content(ligand_id: UUID) -> Optional[LigandContentResponse]:
    return await GetLigandFeature().handle(ligand_id=ligand_id)


@router.post("", summary="Upload ligand")
async def upload_ligand(
        experiment_id: UUID = Form(),
        name: Optional[str] = Form(None),
        smiles: UploadFile = File(None),
        sdf: UploadFile = File(None),
) -> UploadLigandResponse:
    return await UploadLigandFeature().handle(
        request=UploadLigandRequest(
            experiment_id=experiment_id, name=name, smiles=smiles, sdf=sdf
        )
    )


@router.patch("", summary="Update ligand")
async def update_ligand(
        ligand_id: UUID = Form(),
        name: Optional[str] = Form(None),
        smiles: UploadFile = File(None),
        sdf: UploadFile = File(None),
) -> LigandContentResponse:
    return await UpdateLigandFeature().handle(
        request=UpdateLigandRequest(
            ligand_id=ligand_id, name=name, smiles=smiles, sdf=sdf
        )
    )


@router.delete("/{ligand_id}", summary="Delete ligand")
async def delete_ligand(ligand_id: UUID):
    return await DeleteLigandFeature().handle(ligand_id=ligand_id)
