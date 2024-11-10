__all__ = [
    "LigandSearchMetadataQuery",
    "LigandMetadataResponse",
    "LigandSearchContentQuery",
    "LigandContentResponse",
    "UploadLigandRequest",
]

from dataclasses import field
from typing import Optional, Self
from uuid import UUID

from fastapi import UploadFile
from pydantic import BaseModel, model_validator
from pydantic.dataclasses import dataclass


@dataclass
class LigandContentResponse:
    id: UUID
    name: str
    experiment_id: UUID
    smiles_content: Optional[str] = None
    sdf_content: Optional[str] = None
    link: Optional[str] = None
    image: Optional[str] = field(default=None, repr=False)
    drug_likeness: Optional[float] = None
    designed_ligand_score: Optional[float] = None


@dataclass
class LigandMetadataResponse:
    id: UUID
    name: str
    experiment_id: UUID
    smiles_content: Optional[str] = None
    link: Optional[str] = None
    image: Optional[str] = field(default=None, repr=False)
    drug_likeness: Optional[float] = None
    designed_ligand_score: Optional[float] = None


class LigandSearchContentQuery(BaseModel):
    name: Optional[str] = None
    experiment_id: Optional[UUID] = None
    all: Optional[bool] = False

    @model_validator(mode="after")
    def check_at_least_one_field_set(self) -> Self:
        if not self.all and not self.name and not self.experiment_id:
            raise ValueError("You must specify at least one condition to search")

        return self


class LigandSearchMetadataQuery(BaseModel):
    name: Optional[str] = None
    experiment_id: Optional[UUID] = None
    all: Optional[bool] = False

    @model_validator(mode="after")
    def check_at_least_one_field_set(self) -> Self:
        if not self.all and not self.name and not self.experiment_id:
            raise ValueError("You must specify at least one condition to search")

        return self


@dataclass
class UploadLigandRequest:
    experiment_id: UUID
    name: Optional[str] = None
    smiles: Optional[UploadFile] = None
    sdf: Optional[UploadFile] = None
    link: Optional[str] = None


@dataclass
class UpdateLigandRequest:
    ligand_id: UUID
    name: Optional[str] = None
    smiles: Optional[UploadFile] = None
    sdf: Optional[UploadFile] = None
    link: Optional[str] = None


@dataclass
class UploadLigandResponse:
    id: UUID
    name: str
    experiment_id: UUID
    smiles_content: Optional[str] = None
    link: Optional[str] = None
    image: Optional[str] = field(default=None, repr=False)
    drug_likeness: Optional[float] = None
    designed_ligand_score: Optional[float] = None