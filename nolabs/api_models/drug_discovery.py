from pydantic import dataclasses as pcdataclass
import datetime
from fastapi import UploadFile
from typing import List


@pcdataclass.dataclass
class GetExperimentRequest:
    experiment_id: str

@pcdataclass.dataclass
class DeleteExperimentRequest:
    experiment_id: str

@pcdataclass.dataclass
class ExperimentMetadataResponse:
    experiment_id: str
    experiment_name: str
    experiment_date: datetime.datetime

@pcdataclass.dataclass
class ChangeExperimentNameRequest:
    id: str
    name: str

# # #
# Uploading and managing targets section
# # #

@pcdataclass.dataclass
class TargetMetaData:
    target_id: str
    target_name: str

@pcdataclass.dataclass
class UploadTargetRequest:
    experiment_id: str
    fasta_file: UploadFile

@pcdataclass.dataclass
class UploadTargetResponse:
    """
    Returns the list of targets since an uploaded .fasta file could contain multiple sequences inside
    """
    result: List[TargetMetaData]

@pcdataclass.dataclass
class DeleteTargetRequest:
    experiment_id: str
    target_id: str

@pcdataclass.dataclass
class DeleteTargetResponse:
    target_id: str

@pcdataclass.dataclass
class GetTargetDataRequest:
    experiment_id: str
    targe_id: str

@pcdataclass.dataclass
class GetTargetDataResponse:
    protein_sequence: str
    protein_pdb: str | None = None

@pcdataclass.dataclass
class GetTargetsListRequest:
    experiment_id: str

@pcdataclass.dataclass
class GetTargetsListResponse:
    targets: List[TargetMetaData]

@pcdataclass.dataclass
class GetTargetBindingPocketRequest:
    experiment_id: str
    targe_id: str

@pcdataclass.dataclass
class GetTargetBindingPocketResponse:
    pocket_ids: List[int] | None

@pcdataclass.dataclass
class PredictBindingPocketRequest:
    experiment_id: str
    target_id: str
    
@pcdataclass.dataclass
class PredictBindingPocketResponse:
    pocket_ids: List[int] | None

@pcdataclass.dataclass
class PredictMsaRequest:
    experiment_id: str
    target_id: str

@pcdataclass.dataclass
class PredictMsaResponse:
    msa_contents: str

@pcdataclass.dataclass
class GetFoldingRequest:
    experiment_id: str
    targe_id: str

@pcdataclass.dataclass
class GetFoldingResponse:
    pdb_contents: str | None

@pcdataclass.dataclass
class PredictFoldingRequest:
    experiment_id: str
    target_id: str

@pcdataclass.dataclass
class PredictFoldingResponse:
    pdb_content: str | None = None

# # #
# Uploading and managing ligands section
# # #

@pcdataclass.dataclass
class LigandMetaData:
    ligand_id: str
    ligand_name: str

@pcdataclass.dataclass
class UploadLigandRequest:
    experiment_id: str
    target_id: str
    sdf_file: UploadFile

@pcdataclass.dataclass
class UploadLigandResponse:
    ligand_meta_data: LigandMetaData

@pcdataclass.dataclass
class DeleteLigandRequest:
    experiment_id: str
    target_id: str
    ligand_id: str

@pcdataclass.dataclass
class DeleteLigandResponse:
    ligand_id: str

@pcdataclass.dataclass
class GetLigandsListRequest:
    experiment_id: str
    target_id: str

@pcdataclass.dataclass
class GetLigandsListResponse:
    ligands: List[LigandMetaData]

@pcdataclass.dataclass
class GetLigandDataRequest:
    experiment_id: str
    target_id: str
    ligand_id: str

@pcdataclass.dataclass
class GetLigandDataResponse:
    ligand_id: str
    ligand_name: str
    ligand_sdf: str
    ligand_smiles: str

# # #
# Running the experiment and showing results section
# # #
@pcdataclass.dataclass
class ResultMetaData:
    result_id: str
    protein_id: str
    ligand_ids: List[str]

@pcdataclass.dataclass
class DockingRequest:
    experiment_id: str
    target_id: str
    ligand_id: str

@pcdataclass.dataclass
class DockingResponse:
    predicted_pdb: str
    predicted_sdf: str
    plddt_array: List[int]

@pcdataclass.dataclass
class GetResultsListRequest:
    experiment_id: str

@pcdataclass.dataclass
class GetResultsListResponse:
    results_list: List





