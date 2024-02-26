import {
  DrugDiscoveryService,
  TargetMetaData,
  LigandMetaData,
  Body_upload_target_api_v1_drug_discovery_upload_target_post,
  PredictFoldingResponse,
  GetTargetBindingPocketResponse,
  PredictBindingPocketResponse,
  Body_upload_ligand_api_v1_drug_discovery_upload_ligand_post,
  UploadTargetResponse,
  GetTargetDataResponse,
  GetLigandDataResponse,
  ExperimentMetadataResponse,
  OpenAPI,
  type UploadLigandResponse,
  GetTargetMetaDataResponse,
  GetLigandMetaDataResponse,
  RegisterDockingJobResponse,
  RunUmolDockingJobResponse,
  GetUmolDockingResultDataResponse,
  CheckResultDataAvailableResponse,
  CheckJobIsRunningResponse,
  CheckMsaDataAvailableResponse,
  CheckFoldingDataAvailableResponse,
  CheckPocketDataAvailableResponse,
  GetAllResultsListResponse,
  GetResultsListForTargetLigandResponse,
  DeleteDockingJobResponse,
  CheckServiceHealthyResponse,
  GetJobBindingPocketDataResponse,
  GetAllJobsListResponse,
  GetJobsListForTargetLigandResponse,
  RunDiffDockDockingJobResponse,
  GetDiffDockDockingResultDataResponse,
  GetDiffDockLigandSdfResponse,
  GetDiffDockParamsResponse, GetDockingParamsResponse
} from 'src/api/client';
import { CancelablePromise } from 'src/api/client/core/CancelablePromise';
import apiConstants from "src/api/constants";

OpenAPI.BASE = apiConstants.hostname;

// Get a list of experiments
export function getExperimentsApi(): CancelablePromise<Array<ExperimentMetadataResponse>> {
    return DrugDiscoveryService.experimentsApiV1DrugDiscoveryExperimentsGet();
}

// Add a new experiment
export function createExperimentApi(): CancelablePromise<ExperimentMetadataResponse> {
    return DrugDiscoveryService.createExperimentApiV1DrugDiscoveryCreateExperimentGet();
}

// Delete an experiment
export function deleteExperimentApi(experimentId: string): CancelablePromise<any> {
    return DrugDiscoveryService.deleteExperimentApiV1DrugDiscoveryDeleteExperimentDelete(experimentId);
}

// Change the name of an experiment
export function changeExperimentNameApi(experimentId: string, experimentName: string): CancelablePromise<any> {
    return DrugDiscoveryService.changeExperimentNameApiV1DrugDiscoveryChangeExperimentNamePost(experimentId, experimentName);
}

export function getExperimentMetadataApi(experimentId: string): CancelablePromise<any> {
  return DrugDiscoveryService.getExperimentMetadataApiV1DrugDiscoveryGetExperimentMetadataGet(experimentId);
}


// Upload a target
export function uploadTargetApi(formData: Body_upload_target_api_v1_drug_discovery_upload_target_post): CancelablePromise<UploadTargetResponse> {
    return DrugDiscoveryService.uploadTargetApiV1DrugDiscoveryUploadTargetPost(formData);
}

// Delete a target
export function deleteTargetApi(experimentId: string, targetId: string): CancelablePromise<any> {
    return DrugDiscoveryService.deleteTargetApiV1DrugDiscoveryDeleteTargetDelete(experimentId, targetId);
}

// Get targets list
export function getTargetsListApi(experimentId: string): CancelablePromise<Array<TargetMetaData>> {
    return DrugDiscoveryService.getTargetsListApiV1DrugDiscoveryGetTargetsListGet(experimentId);
}

export function getTargetMetaDataApi(experimentId: string, targetId: string): CancelablePromise<GetTargetMetaDataResponse> {
  return DrugDiscoveryService.getTargetMetaDataApiV1DrugDiscoveryGetTargetMetaDataGet(experimentId, targetId);
}

export function changeTargetNameApi(experimentId: string, targetId: string, targetName: string): CancelablePromise<GetTargetMetaDataResponse> {
  return DrugDiscoveryService.updateTargetNameApiV1DrugDiscoveryUpdateTargetNamePost(experimentId, targetId, targetName);
}

export function getTargetDataApi(experimentId: string, targetId: string): CancelablePromise<GetTargetDataResponse> {
    return DrugDiscoveryService.getTargetDataApiV1DrugDiscoveryGetTargetDataGet(experimentId, targetId);
}

// Upload a ligand
export function uploadLigandApi(formData: Body_upload_ligand_api_v1_drug_discovery_upload_ligand_post): CancelablePromise<UploadLigandResponse> {
    return DrugDiscoveryService.uploadLigandApiV1DrugDiscoveryUploadLigandPost(formData);
}

// Delete a ligand
export function deleteLigandApi(experimentId: string, targetId: string, ligandId: string): CancelablePromise<any> {
    return DrugDiscoveryService.deleteLigandApiV1DrugDiscoveryDeleteLigandDelete(experimentId, targetId, ligandId);
}

// Get ligands list
export function getLigandsListApi(experimentId: string, targetId: string): CancelablePromise<Array<LigandMetaData>> {
    return DrugDiscoveryService.getLigandsListApiV1DrugDiscoveryGetLigandsListGet(experimentId, targetId);
}

export function getLigandMetaDataApi(experimentId: string, targetId: string, ligandId: string): CancelablePromise<GetLigandMetaDataResponse> {
  return DrugDiscoveryService.getLigandDataApiV1DrugDiscoveryGetLigandMetaDataGet(experimentId, targetId, ligandId);
}

export function getLigandDataApi(experimentId: string, targetId: string, ligandId: string): CancelablePromise<GetLigandDataResponse> {
    return DrugDiscoveryService.getLigandDataApiV1DrugDiscoveryGetLigandDataGet(experimentId, targetId, ligandId);
}

// Get target binding pocket
export function getTargetBindingPocketApi(experimentId: string, targetId: string, foldingMethod: string): CancelablePromise<GetTargetBindingPocketResponse> {
    return DrugDiscoveryService.getTargetBindingPocketApiV1DrugDiscoveryGetTargetBindingPocketGet(experimentId, targetId, foldingMethod);
}

export function setTargetBindingPocketApi(experimentId: string, targetId: string, pocketIds: Array<number>) {
  return DrugDiscoveryService.setTargetBindingPocketApiV1DrugDiscoverySetTargetBindingPocketPost(experimentId, targetId, pocketIds);
}


// Predict binding pocket
export function predictBindingPocketApi(experimentId: string, targetId: string, foldingMethod: string): CancelablePromise<PredictBindingPocketResponse> {
    return DrugDiscoveryService.predictBindingPocketApiV1DrugDiscoveryPredictBindingPocketPost(experimentId, targetId, foldingMethod);
}

// Predict folding
export function predictLightFoldingApi(experimentId: string, targetId: string): CancelablePromise<PredictFoldingResponse> {
    return DrugDiscoveryService.predictFoldingApiV1DrugDiscoveryPredictEsmfoldLightPost(experimentId, targetId);
}

export function predictFoldingApi(experimentId: string, targetId: string): CancelablePromise<PredictFoldingResponse> {
  return DrugDiscoveryService.predictFoldingApiV1DrugDiscoveryPredictEsmfoldPost(experimentId, targetId);
}

export function registerDockingJobApi(experimentId: string, targetId: string, ligandId: string, foldingMethod: string): CancelablePromise<RegisterDockingJobResponse> {
  return DrugDiscoveryService.registerDockingJobApiV1DrugDiscoveryRegisterDockingJobPost(experimentId, targetId, ligandId, foldingMethod);
}

export function runUmolDockingJobApi(experimentId: string, targetId: string, ligandId: string, jobId: string): CancelablePromise<RunUmolDockingJobResponse> {
  return DrugDiscoveryService.runUmolDockingApiV1DrugDiscoveryRunUmolDockingJobPost(experimentId, targetId, ligandId, jobId);
}

export function runDiffDockDockingJobApi(experimentId: string, targetId: string, ligandId: string, jobId: string): CancelablePromise<RunDiffDockDockingJobResponse> {
  return DrugDiscoveryService.runDiffdockDockingApiV1DrugDiscoveryRunDiffdockDockingJobPost(experimentId, targetId, ligandId, jobId);
}


export function getUmolDockingJobResultDataApi(experimentId: string, targetId: string, ligandId: string, jobId: string): CancelablePromise<GetUmolDockingResultDataResponse> {
  return DrugDiscoveryService.getUmolDockingResultDataApiV1DrugDiscoveryGetUmolDockingResultDataGet(experimentId, targetId, ligandId, jobId);
}

export function getDiffDockParamsApi(experimentId: string, targetId: string, ligandId: string, jobId: string): CancelablePromise<GetDiffDockParamsResponse> {
  return DrugDiscoveryService.getDiffdockParamsApiV1DrugDiscoveryGetDiffdockParamsGet(experimentId, targetId, ligandId, jobId);
}

export function updateDockingParamsApi(experimentId: string, targetId: string, ligandId: string, jobId: string, foldingMethod: string, dockingMethod: string): CancelablePromise<GetDockingParamsResponse> {
  return DrugDiscoveryService.updateDockingParamsApiV1DrugDiscoveryUpdateDockingParamsPost(experimentId, targetId, ligandId, jobId, foldingMethod, dockingMethod);
}


export function updateDiffDockParamsApi(experimentId: string, targetId: string, ligandId: string, jobId: string, samples_per_complex: number): CancelablePromise<GetDiffDockParamsResponse> {
  return DrugDiscoveryService.updateDiffdockParamsApiV1DrugDiscoveryUpdateDiffdockParamsPost(experimentId, targetId, ligandId, jobId, samples_per_complex);
}

export function getDiffDockDockingJobResultDataApi(experimentId: string, targetId: string, ligandId: string, jobId: string): CancelablePromise<GetDiffDockDockingResultDataResponse> {
  return DrugDiscoveryService.getDiffdockDockingResultDataApiV1DrugDiscoveryGetDiffdockDockingResultDataGet(experimentId, targetId, ligandId, jobId);
}

export function getDiffDockLigandSdfApi(experimentId: string, targetId: string, ligandId: string, jobId: string, ligandFileName: string): CancelablePromise<GetDiffDockLigandSdfResponse> {
  return DrugDiscoveryService.getDiffdockLigandSdfApiV1DrugDiscoveryGetDiffdockLigandSdfGet(experimentId, targetId, ligandId, jobId, ligandFileName);
}


export function checkDockingResultAvailableApi(experimentId: string, targetId: string, ligandId: string, jobId: string): CancelablePromise<CheckResultDataAvailableResponse> {
  return DrugDiscoveryService.checkResultDataAvailableApiV1DrugDiscoveryCheckResultDataAvailableGet(experimentId, targetId, ligandId, jobId);
}

export function checkMsaDataAvailableApi(experimentId: string, targetId: string): CancelablePromise<CheckMsaDataAvailableResponse> {
  return DrugDiscoveryService.checkMsaDataAvailableApiV1DrugDiscoveryCheckMsaDataAvailableGet(experimentId, targetId);
}

export function checkMsaJobIsRunningApi(jobId: string): CancelablePromise<CheckJobIsRunningResponse> {
  return DrugDiscoveryService.checkMsaJobRunningApiV1DrugDiscoveryCheckMsaJobRunningGet(jobId);
}

export function checkFoldingDataAvailableApi(experimentId: string, targetId: string, foldingMethod: string): CancelablePromise<CheckFoldingDataAvailableResponse> {
  return DrugDiscoveryService.checkFoldingDataAvailableApiV1DrugDiscoveryCheckFoldingDataAvailableGet(experimentId, targetId, foldingMethod);
}

export function checkEsmFoldJobIsRunningApi(jobId: string): CancelablePromise<CheckJobIsRunningResponse> {
  return DrugDiscoveryService.checkEsmfoldJobRunningApiV1DrugDiscoveryCheckEsmfoldJobRunningGet(jobId);
}

export function checkEsmFoldLightJobIsRunningApi(jobId: string): CancelablePromise<CheckJobIsRunningResponse> {
  return DrugDiscoveryService.checkEsmfoldLightJobRunningApiV1DrugDiscoveryCheckEsmfoldLightJobRunningGet(jobId);
}

export function checkPocketDataAvailableApi(experimentId: string, targetId: string, ligandId: string, jobId: string): CancelablePromise<CheckPocketDataAvailableResponse> {
  return DrugDiscoveryService.checkPocketDataAvailableApiV1DrugDiscoveryCheckPocketDataAvailableGet(experimentId, targetId, ligandId, jobId);
}

export function checkP2RankJobIsRunningApi(jobId: string): CancelablePromise<CheckJobIsRunningResponse> {
  return DrugDiscoveryService.checkP2RankJobRunningApiV1DrugDiscoveryCheckP2RankJobRunningGet(jobId);
}

export function checkUmolJobIsRunningApi(jobId: string): CancelablePromise<CheckJobIsRunningResponse> {
  return DrugDiscoveryService.checkUmolJobRunningApiV1DrugDiscoveryCheckUmolJobRunningGet(jobId);
}

export function checkDiffDockJobIsRunningApi(jobId: string): CancelablePromise<CheckJobIsRunningResponse> {
  return DrugDiscoveryService.checkDiffdockJobRunningApiV1DrugDiscoveryCheckDiffdockJobRunningGet(jobId);
}

export function getAllDockingResultsListApi(experimentId: string): CancelablePromise<GetAllResultsListResponse> {
  return DrugDiscoveryService.getAllResultsListApiV1DrugDiscoveryGetAllResultsListGet(experimentId);
}

export function getAllDockingJobsListApi(experimentId: string): CancelablePromise<GetAllJobsListResponse> {
  return DrugDiscoveryService.getAllJobsListApiV1DrugDiscoveryGetAllJobsListGet(experimentId);
}

export function getDockingResultsListForTargetLigandApi(experimentId: string, targetId: string, ligandId: string): CancelablePromise<GetResultsListForTargetLigandResponse> {
  return DrugDiscoveryService.getResultsListForTargetLigandApiV1DrugDiscoveryGetResultsListForTargetLigandGet(experimentId, targetId, ligandId);
}

export function getDockingJobsListForTargetLigandApi(experimentId: string, targetId: string, ligandId: string): CancelablePromise<GetJobsListForTargetLigandResponse> {
  return DrugDiscoveryService.getJobsListForTargetLigandApiV1DrugDiscoveryGetJobsListForTargetLigandGet(experimentId, targetId, ligandId);
}

export function deleteDockingJobApi(experimentId: string, targetId: string, ligandId: string, jobId: string): CancelablePromise<DeleteDockingJobResponse> {
  return DrugDiscoveryService.deleteDockingJobApiV1DrugDiscoveryDeleteDockingJobDelete(experimentId, targetId, ligandId, jobId);
}

export function checkMsaServiceHealthApi(): CancelablePromise<CheckServiceHealthyResponse> {
  return DrugDiscoveryService.checkMsaServiceHealthApiV1DrugDiscoveryCheckMsaServiceHealthGet();
}

export function checkP2RankServiceHealthApi(): CancelablePromise<CheckServiceHealthyResponse> {
  return DrugDiscoveryService.checkP2RankServiceHealthApiV1DrugDiscoveryCheckP2RankServiceHealthGet();
}

export function checkEsmFoldServiceHealthApi(): CancelablePromise<CheckServiceHealthyResponse> {
  return DrugDiscoveryService.checkEsmfoldServiceHealthApiV1DrugDiscoveryCheckEsmfoldServiceHealthGet();
}

export function checkEsmFoldLightServiceHealthApi(): CancelablePromise<CheckServiceHealthyResponse> {
  return DrugDiscoveryService.checkEsmfoldLightServiceHealthApiV1DrugDiscoveryCheckEsmfoldLightServiceHealthGet();
}

export function checkUmolServiceHealthApi(): CancelablePromise<CheckServiceHealthyResponse> {
  return DrugDiscoveryService.checkUmolServiceHealthApiV1DrugDiscoveryCheckUmolServiceHealthGet();
}

export function checkDiffDockServiceHealthApi(): CancelablePromise<CheckServiceHealthyResponse> {
  return DrugDiscoveryService.checkDiffdockServiceHealthApiV1DrugDiscoveryCheckDiffdockServiceHealthGet();
}

export function getJobPocketDataApi(experimentId: string, targetId: string, ligandId: string, jobId: string): CancelablePromise<GetJobBindingPocketDataResponse> {
  return DrugDiscoveryService.getJobPocketDataApiV1DrugDiscoveryGetJobPocketDataGet(experimentId, targetId, ligandId, jobId);
}
