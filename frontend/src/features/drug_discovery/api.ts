import {
    DrugDiscoveryService,
    nolabs__api_models__drug_discovery__ExperimentMetadataResponse,
    TargetMetaData,
    LigandMetaData,
    ChangeExperimentNameRequest,
    Body_upload_target_api_v1_drug_discovery_upload_target_post,
    DockingRequest,
    GetFoldingRequest,
    PredictBindingPocketRequest,
    PredictFoldingRequest,
    PredictMsaRequest,
    Body_upload_ligand_api_v1_drug_discovery_upload_ligand_post
} from '../../api/client';
import { CancelablePromise } from 'api/client/core/CancelablePromise';

// Get a list of experiments
export function getExperiments(): CancelablePromise<Array<nolabs__api_models__drug_discovery__ExperimentMetadataResponse>> {
    return DrugDiscoveryService.experimentsApiV1DrugDiscoveryExperimentsGet();
}

// Add a new experiment
export function addExperiment(): CancelablePromise<nolabs__api_models__drug_discovery__ExperimentMetadataResponse> {
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
export function uploadTarget(formData: Body_upload_target_api_v1_drug_discovery_upload_target_post): CancelablePromise<any> {
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

// Get target binding pocket
export function getTargetBindingPocket(experimentId: string, targetId: string): CancelablePromise<any> {
    return DrugDiscoveryService.getTargetBindingPocketApiV1DrugDiscoveryGetTargetBindingPocketGet(experimentId, targetId);
}

// Predict binding pocket
export function predictBindingPocket(requestBody: PredictBindingPocketRequest): CancelablePromise<any> {
    return DrugDiscoveryService.predictBindingPocketApiV1DrugDiscoveryPredictBindingPocketPost(requestBody);
}

// Predict MSA (Multiple Sequence Alignment)
export function predictMsa(requestBody: PredictMsaRequest): CancelablePromise<any> {
    return DrugDiscoveryService.predictMsaApiV1DrugDiscoveryPredictMsaPost(requestBody);
}

// Check if a structure has been folded
export function checkFoldingExist(requestBody: GetFoldingRequest): CancelablePromise<any> {
    return DrugDiscoveryService.getFoldedStructureApiV1DrugDiscoveryGetFoldedStructureGet(requestBody);
}

// Predict folding
export function predictFolding(requestBody: PredictFoldingRequest): CancelablePromise<any> {
    return DrugDiscoveryService.predictFoldingApiV1DrugDiscoveryPredictFoldingPost(requestBody);
}

// Perform docking
export function performDocking(requestBody: DockingRequest): CancelablePromise<any> {
    return DrugDiscoveryService.performDockingApiV1DrugDiscoveryPredictDockingPost(requestBody);
}