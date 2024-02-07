from pydantic import dataclasses as pcdataclass
from fastapi import UploadFile
from typing import List


@pcdataclass.dataclass
class GetExperimentRequest:
    experiment_id: str

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
class GetTargetMetaDataRequest:
    experiment_id: str
    target_id: str

@pcdataclass.dataclass
class GetTargetMetaDataResponse:
    target_id: str
    target_name: str

@pcdataclass.dataclass
class GetTargetDataRequest:
    experiment_id: str
    target_id: str

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
    target_id: str

@pcdataclass.dataclass
class GetTargetBindingPocketResponse:
    pocket_ids: List[int] | None

@pcdataclass.dataclass
class SetTargetBindingPocketRequest:
    experiment_id: str
    target_id: str
    pocket_ids: List[int]

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
class CheckPocketDataAvailableRequest:
    experiment_id: str
    target_id: str
    ligand_id: str
    job_id: str


@pcdataclass.dataclass
class CheckPocketDataAvailableResponse:
    is_available: bool

@pcdataclass.dataclass
class GetJobBindingPocketDataRequest:
    experiment_id: str
    target_id: str
    ligand_id: str
    job_id: str

@pcdataclass.dataclass
class GetJobBindingPocketDataResponse:
    pocket_ids: List[int] | None


@pcdataclass.dataclass
class CheckMsaDataAvailableRequest:
    experiment_id: str
    target_id: str


@pcdataclass.dataclass
class CheckMsaDataAvailableResponse:
    is_available: bool


@pcdataclass.dataclass
class CheckFoldingDataAvailableRequest:
    experiment_id: str
    target_id: str


@pcdataclass.dataclass
class CheckFoldingDataAvailableResponse:
    is_available: bool

@pcdataclass.dataclass
class GetFoldingRequest:
    experiment_id: str
    target_id: str

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
class GetLigandMetaDataRequest:
    experiment_id: str
    target_id: str
    ligand_id: str

@pcdataclass.dataclass
class GetLigandMetaDataResponse:
    ligand_id: str
    ligand_name: str

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
    experiment_id: str
    job_id: str
    target_id: str
    ligand_id: str

@pcdataclass.dataclass
class RegisterDockingJobRequest:
    experiment_id: str
    target_id: str
    ligand_id: str

@pcdataclass.dataclass
class RegisterDockingJobResponse:
    job_id: str

@pcdataclass.dataclass
class DeleteDockingJobRequest:
    experiment_id: str
    target_id: str
    ligand_id: str
    job_id: str

@pcdataclass.dataclass
class DeleteDockingJobResponse:
    job_id: str

@pcdataclass.dataclass
class RunDockingJobRequest:
    experiment_id: str
    target_id: str
    ligand_id: str
    job_id: str

@pcdataclass.dataclass
class RunDockingJobResponse:
    predicted_pdb: str
    predicted_sdf: str
    plddt_array: List[int]
    job_id: str


@pcdataclass.dataclass
class CheckResultDataAvailableRequest:
    experiment_id: str
    target_id: str
    ligand_id: str
    job_id: str

@pcdataclass.dataclass
class CheckResultDataAvailableResponse:
    result_available: bool

@pcdataclass.dataclass
class GetDockingResultDataRequest:
    experiment_id: str
    target_id: str
    ligand_id: str
    job_id: str

@pcdataclass.dataclass
class GetDockingResultDataResponse:
    predicted_pdb: str
    predicted_sdf: str
    plddt_array: List[int]
    job_id: str

@pcdataclass.dataclass
class CheckJobIsRunningRequest:
    job_id: str

@pcdataclass.dataclass
class CheckJobIsRunningResponse:
    is_running: bool

@pcdataclass.dataclass
class CheckResultAvailableRequest:
    experiment_id: str
    target_id: str
    ligandId: str
    job_id: str

@pcdataclass.dataclass
class GetAllResultsListRequest:
    experiment_id: str

@pcdataclass.dataclass
class GetAllResultsListResponse:
    results_list: List[ResultMetaData]

@pcdataclass.dataclass
class GetResultsListForTargetLigandRequest:
    experiment_id: str
    target_id: str
    ligand_id: str

@pcdataclass.dataclass
class GetResultsListForTargetLigandResponse:
    results_list: List[ResultMetaData]

@pcdataclass.dataclass
class CheckServiceHealthyResponse:
    is_healthy: bool