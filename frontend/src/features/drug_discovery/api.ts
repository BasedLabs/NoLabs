import {
    DrugDiscoveryService,
    TargetMetaData,
    LigandMetaData,
    ChangeExperimentNameRequest,
    Body_upload_target_api_v1_drug_discovery_upload_target_post,
    DockingRequest,
    GetFoldingRequest,
    PredictFoldingResponse,
    PredictBindingPocketResponse,
    Body_upload_ligand_api_v1_drug_discovery_upload_ligand_post,
    UploadTargetResponse,
    GetTargetDataResponse,
    GetLigandDataResponse, ExperimentMetadataResponse, OpenAPI
} from 'src/api/client';
import { CancelablePromise } from 'src/api/client/core/CancelablePromise';
import apiConstants from "src/api/constants";

OpenAPI.BASE = apiConstants.hostname;

// Get a list of experiments
export function getExperiments(): CancelablePromise<Array<ExperimentMetadataResponse>> {
    return DrugDiscoveryService.experimentsApiV1DrugDiscoveryExperimentsGet();
}

// Add a new experiment
export function addExperiment(): CancelablePromise<ExperimentMetadataResponse> {
    return DrugDiscoveryService.addExperimentApiV1DrugDiscoveryAddExperimentPost();
}

// Delete an experiment
export function deleteExperiment(experimentId: string): CancelablePromise<any> {
    return DrugDiscoveryService.deleteExperimentApiV1DrugDiscoveryDeleteExperimentDelete(experimentId);
}

// Change the name of an experiment
export function changeExperimentName(requestBody: ChangeExperimentNameRequest): CancelablePromise<any> {
    return DrugDiscoveryService.changeExperimentNameApiV1DrugDiscoveryChangeExperimentNamePost(requestBody);
}

// Upload a target
export function uploadTarget(formData: Body_upload_target_api_v1_drug_discovery_upload_target_post): CancelablePromise<UploadTargetResponse> {
    return DrugDiscoveryService.uploadTargetApiV1DrugDiscoveryUploadTargetPost(formData);
}

// Delete a target
export function deleteTarget(experimentId: string, targetId: string): CancelablePromise<any> {
    return DrugDiscoveryService.deleteTargetApiV1DrugDiscoveryDeleteTargetDelete(experimentId, targetId);
}

// Get targets list
export function getTargetsList(experimentId: string): CancelablePromise<Array<TargetMetaData>> {
    return DrugDiscoveryService.getTargetsListApiV1DrugDiscoveryGetTargetsListGet(experimentId);
}

export function getTargetData(experimentId: string, targetId: string): CancelablePromise<GetTargetDataResponse> {
    return DrugDiscoveryService.getTargetsListApiV1DrugDiscoveryGetTargetDataGet(experimentId, targetId);
}

// Upload a ligand
export function uploadLigand(formData: Body_upload_ligand_api_v1_drug_discovery_upload_ligand_post): CancelablePromise<any> {
    return DrugDiscoveryService.uploadLigandApiV1DrugDiscoveryUploadLigandPost(formData);
}

// Delete a ligand
export function deleteLigand(experimentId: string, targetId: string, ligandId: string): CancelablePromise<any> {
    return DrugDiscoveryService.deleteLigandApiV1DrugDiscoveryDeleteLigandDelete(experimentId, targetId, ligandId);
}

// Get ligands list
export function getLigandsList(experimentId: string, targetId: string): CancelablePromise<Array<LigandMetaData>> {
    return DrugDiscoveryService.getLigandsListApiV1DrugDiscoveryGetLigandsListGet(experimentId, targetId);
}

export function getLigandData(experimentId: string, targetId: string, ligandId: string): CancelablePromise<GetLigandDataResponse> {
    return DrugDiscoveryService.getLigandDataApiV1DrugDiscoveryGetLigandDataGet(experimentId, targetId, ligandId);
}

// Get target binding pocket
export function getTargetBindingPocket(experimentId: string, targetId: string): CancelablePromise<any> {
    return DrugDiscoveryService.getTargetBindingPocketApiV1DrugDiscoveryGetTargetBindingPocketGet(experimentId, targetId);
}

// Predict binding pocket
export function predictBindingPocket(experimentId: string, targetId: string): CancelablePromise<PredictBindingPocketResponse> {
    return DrugDiscoveryService.predictBindingPocketApiV1DrugDiscoveryPredictBindingPocketPost(experimentId, targetId);
}

// Predict MSA (Multiple Sequence Alignment)
export function predictMsa(experimentId: string, targetId: string): CancelablePromise<any> {
    return DrugDiscoveryService.predictMsaApiV1DrugDiscoveryPredictMsaPost(experimentId, targetId);
}

// Check if a structure has been folded
export function checkFoldingExist(requestBody: GetFoldingRequest): CancelablePromise<any> {
    return DrugDiscoveryService.getFoldedStructureApiV1DrugDiscoveryGetFoldedStructureGet(requestBody);
}

// Predict folding
export function predictFolding(experimentId: string, targetId: string): CancelablePromise<PredictFoldingResponse> {
    return DrugDiscoveryService.predictFoldingApiV1DrugDiscoveryPredictFoldingPost(experimentId, targetId);
}

// Perform docking
export function performDocking(requestBody: DockingRequest): CancelablePromise<any> {
    return DrugDiscoveryService.performDockingApiV1DrugDiscoveryPredictDockingPost(requestBody);
}