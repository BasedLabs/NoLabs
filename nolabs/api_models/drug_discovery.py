from pydantic import dataclasses as pcdataclass
from fastapi import UploadFile
from typing import List, Optional


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
    folding_method: Optional[str] = None

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
class UpdateTargetNameRequest:
    experiment_id: str
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
    folding_method: str


@pcdataclass.dataclass
class CheckFoldingDataAvailableResponse:
    is_available: bool

@pcdataclass.dataclass
class GetFoldingRequest:
    experiment_id: str
    target_id: str
    folding_method: str

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
class UploadTargetLigandRequest:
    experiment_id: str
    target_id: str
    sdf_file: UploadFile

@pcdataclass.dataclass
class UploadTargetLigandResponse:
    ligand_meta_data: LigandMetaData

@pcdataclass.dataclass
class UploadLoneLigandRequest:
    experiment_id: str
    sdf_file: UploadFile

@pcdataclass.dataclass
class UploadLoneLigandResponse:
    ligand_meta_data: LigandMetaData

@pcdataclass.dataclass
class DeleteTargetLigandRequest:
    experiment_id: str
    target_id: str
    ligand_id: str

@pcdataclass.dataclass
class DeleteTargetLigandResponse:
    ligand_id: str

@pcdataclass.dataclass
class DeleteLoneLigandRequest:
    experiment_id: str
    ligand_id: str

@pcdataclass.dataclass
class DeleteLoneLigandResponse:
    ligand_id: str

@pcdataclass.dataclass
class GetTargetLigandsListRequest:
    experiment_id: str
    target_id: str

@pcdataclass.dataclass
class GetTargetLigandsListResponse:
    ligands: List[LigandMetaData]

@pcdataclass.dataclass
class GetLoneLigandsListRequest:
    experiment_id: str

@pcdataclass.dataclass
class GetLoneLigandsListResponse:
    ligands: List[LigandMetaData]

@pcdataclass.dataclass
class GetTargetLigandMetaDataRequest:
    experiment_id: str
    target_id: str
    ligand_id: str

@pcdataclass.dataclass
class GetTargetLigandMetaDataResponse:
    ligand_id: str
    ligand_name: str

@pcdataclass.dataclass
class GetLoneLigandMetaDataRequest:
    experiment_id: str
    ligand_id: str

@pcdataclass.dataclass
class GetLoneLigandMetaDataResponse:
    ligand_id: str
    ligand_name: str

@pcdataclass.dataclass
class GetTargetLigandDataRequest:
    experiment_id: str
    target_id: str
    ligand_id: str

@pcdataclass.dataclass
class GetTargetLigandDataResponse:
    ligand_id: str
    ligand_name: str
    ligand_sdf: str
    ligand_smiles: str

@pcdataclass.dataclass
class GetLoneLigandDataRequest:
    experiment_id: str
    ligand_id: str

@pcdataclass.dataclass
class GetLoneLigandDataResponse:
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
class JobMetaData:
    experiment_id: str
    job_id: str
    target_id: str
    ligand_id: str
    folding_method: str
    docking_method: str

@pcdataclass.dataclass
class RegisterDockingJobRequest:
    experiment_id: str
    target_id: str
    ligand_id: str
    folding_method: str

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
class RunUmolDockingJobRequest:
    experiment_id: str
    target_id: str
    ligand_id: str
    job_id: str

@pcdataclass.dataclass
class RunUmolDockingJobResponse:
    predicted_pdb: str
    predicted_sdf: str
    plddt_array: List[int]
    job_id: str

@pcdataclass.dataclass
class RunDiffDockDockingJobRequest:
    experiment_id: str
    target_id: str
    ligand_id: str
    job_id: str

@pcdataclass.dataclass
class DiffDockLigandMetaData:
    job_id: str
    target_id: str
    ligand_id: str
    predicted_ligand_file_name: str
    minimized_affinity: float
    scored_affinity: float
    confidence: float | None = None

@pcdataclass.dataclass
class RunDiffDockDockingJobResponse:
    predicted_pdb: str
    predicted_ligands: List[DiffDockLigandMetaData]

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
class GetUmolDockingResultDataResponse:
    predicted_pdb: str
    predicted_sdf: str
    plddt_array: List[int]
    job_id: str

@pcdataclass.dataclass
class GetDockingParamsRequest:
    experiment_id: str
    target_id: str
    ligand_id: str
    job_id: str

@pcdataclass.dataclass
class GetDockingParamsResponse:
    folding_method: str
    docking_method: str

@pcdataclass.dataclass
class UpdateDockingParamsRequest:
    experiment_id: str
    target_id: str
    ligand_id: str
    job_id: str
    folfing_method: str
    docking_method: str



@pcdataclass.dataclass
class GetDiffDockParamsRequest:
    experiment_id: str
    target_id: str
    ligand_id: str
    job_id: str

@pcdataclass.dataclass
class GetDiffDockParamsResponse:
    samples_per_complex: int


@pcdataclass.dataclass
class UpdateDiffDockParamsRequest:
    experiment_id: str
    target_id: str
    ligand_id: str
    job_id: str
    samples_per_complex: int


@pcdataclass.dataclass
class GetDiffDockDockingResultDataResponse:
    predicted_pdb: str
    predicted_ligands: List[DiffDockLigandMetaData]

@pcdataclass.dataclass
class GetDiffDockLigandSdfRequest:
    experiment_id: str
    target_id: str
    ligand_id: str
    job_id: str
    ligand_file_name: str

@pcdataclass.dataclass
class GetDiffDockLigandSdfResponse:
    sdf_contents: str


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
    results_list: List[JobMetaData]

@pcdataclass.dataclass
class GetAllJobsListRequest:
    experiment_id: str

@pcdataclass.dataclass
class GetAllJobsListResponse:
    jobs_list: List[JobMetaData]

@pcdataclass.dataclass
class GetResultsListForTargetLigandRequest:
    experiment_id: str
    target_id: str
    ligand_id: str

@pcdataclass.dataclass
class GetResultsListForTargetLigandResponse:
    results_list: List[JobMetaData]

@pcdataclass.dataclass
class GetJobsListForTargetLigandRequest:
    experiment_id: str
    target_id: str
    ligand_id: str

@pcdataclass.dataclass
class GetJobsListForTargetLigandResponse:
    jobs_list: List[JobMetaData]

@pcdataclass.dataclass
class CheckServiceHealthyResponse:
    is_healthy: bool