from pydantic import dataclasses as pcdataclass
import datetime
from fastapi import UploadFile
from typing import List, Optional, Dict


@pcdataclass.dataclass
class GetExperimentRequest:
    experimentId: str

@pcdataclass.dataclass
class DeleteExperimentRequest:
    id: str

@pcdataclass.dataclass
class ExperimentMetadataResponse:
    id: str
    name: str
    date: datetime.datetime

@pcdataclass.dataclass
class ChangeExperimentNameRequest:
    id: str
    name: str

# # #
# Uploading and managing targets section
# # #

@pcdataclass.dataclass
class TargetMetaData:
    targetId: str
    targetName: str

@pcdataclass.dataclass
class UploadTargetRequest:
    experimentId: str
    fasta_file: UploadFile

@pcdataclass.dataclass
class UploadTargetResponse:
    """
    Returns the list of targets since an uploaded .fasta file could contain multiple sequences inside
    """
    result: List[TargetMetaData]

@pcdataclass.dataclass
class GetTargetDataRequest:
    experimentId: str
    targetId: str

@pcdataclass.dataclass
class GetTargetDataResponse:
    protein_pdb: str | None = None
    protein_fasta: str

@pcdataclass.dataclass
class GetTargetsListRequest:
    experimentId: str

@pcdataclass.dataclass
class GetTargetsListResponse:
    targets: List[UploadTargetResponse]

@pcdataclass.dataclass
class GetTargetBindingPocketRequest:
    experimentId: str
    targetId: str

@pcdataclass.dataclass
class GetTargetBindingPocketResponse:
    pocketIds: List[int] | None

@pcdataclass.dataclass
class PredictBindingPocketRequest:
    experimentId: str
    targetId: str
    
@pcdataclass.dataclass
class PredictBindingPocketResponse:
    pocketIds: List[int] | None

# # #
# Uploading and managing ligands section
# # #

@pcdataclass.dataclass
class LigandMetaData:
    ligandId: str
    ligandName: str

@pcdataclass.dataclass
class UploadLigandRequest:
    experimentId: str
    sdf: UploadFile

@pcdataclass.dataclass
class UploadLigandResponse:
    ligandId: str
    ligandMetaData: LigandMetaData

@pcdataclass.dataclass
class GetLigandsListRequest:
    experimentId: str

@pcdataclass.dataclass
class GetLigandsListResponse:
    targets: List[UploadLigandResponse]

@pcdataclass.dataclass
class GetLigandDataRequest:
    ligandId: str

@pcdataclass.dataclass
class GetLigandDataResponse:
    ligand_sdf: str

# # #
# Running the experiment and showing results section
# # #
@pcdataclass.dataclass
class ResultMetaData:
    resultId: str
    proteinId: str
    ligandIds: List[str]

@pcdataclass.dataclass
class DockingRequest:
    experimentId: str
    protein_id: str
    ligand_id: str

@pcdataclass.dataclass
class DockingResponse:
    predicted_pdb: str
    predicted_sdf: str

@pcdataclass.dataclass
class GetResultsListRequest:
    experimentId: str

@pcdataclass.dataclass
class GetResultsListResponse:
    resultsList: List





