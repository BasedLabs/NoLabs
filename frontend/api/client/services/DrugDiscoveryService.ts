/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { Body_upload_ligand_api_v1_drug_discovery_upload_ligand_post } from '../models/Body_upload_ligand_api_v1_drug_discovery_upload_ligand_post';
import type { Body_upload_target_api_v1_drug_discovery_upload_target_post } from '../models/Body_upload_target_api_v1_drug_discovery_upload_target_post';
import type { ChangeExperimentNameRequest } from '../models/ChangeExperimentNameRequest';
import type { DockingRequest } from '../models/DockingRequest';
import type { GetFoldingRequest } from '../models/GetFoldingRequest';
import type { nolabs__api_models__drug_discovery__ExperimentMetadataResponse } from '../models/nolabs__api_models__drug_discovery__ExperimentMetadataResponse';
import type { PredictBindingPocketRequest } from '../models/PredictBindingPocketRequest';
import type { PredictFoldingRequest } from '../models/PredictFoldingRequest';
import type { PredictMsaRequest } from '../models/PredictMsaRequest';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class DrugDiscoveryService {
    /**
     * Experiments
     * @returns nolabs__api_models__drug_discovery__ExperimentMetadataResponse Successful Response
     * @throws ApiError
     */
    public static experimentsApiV1DrugDiscoveryExperimentsGet(): CancelablePromise<Array<nolabs__api_models__drug_discovery__ExperimentMetadataResponse>> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/experiments',
        });
    }
    /**
     * Add Experiment
     * @returns nolabs__api_models__drug_discovery__ExperimentMetadataResponse Successful Response
     * @throws ApiError
     */
    public static addExperimentApiV1DrugDiscoveryAddExperimentPost(): CancelablePromise<nolabs__api_models__drug_discovery__ExperimentMetadataResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/add-experiment',
        });
    }
    /**
     * Delete Experiment
     * @param experimentId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static deleteExperimentApiV1DrugDiscoveryDeleteExperimentDelete(
        experimentId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'DELETE',
            url: '/api/v1/drug_discovery/delete-experiment',
            query: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Change Experiment Name
     * @param requestBody
     * @returns any Successful Response
     * @throws ApiError
     */
    public static changeExperimentNameApiV1DrugDiscoveryChangeExperimentNamePost(
        requestBody: ChangeExperimentNameRequest,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/change-experiment-name',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Upload Target
     * @param formData
     * @returns any Successful Response
     * @throws ApiError
     */
    public static uploadTargetApiV1DrugDiscoveryUploadTargetPost(
        formData: Body_upload_target_api_v1_drug_discovery_upload_target_post,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/upload-target',
            formData: formData,
            mediaType: 'multipart/form-data',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Delete Target
     * @param experimentId
     * @param targetId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static deleteTargetApiV1DrugDiscoveryDeleteTargetDelete(
        experimentId: string,
        targetId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'DELETE',
            url: '/api/v1/drug_discovery/delete-target',
            query: {
                'experiment_id': experimentId,
                'target_id': targetId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get Targets List
     * @param experimentId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static getTargetsListApiV1DrugDiscoveryGetTargetsListGet(
        experimentId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/get-targets-list',
            query: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Upload Ligand
     * @param formData
     * @returns any Successful Response
     * @throws ApiError
     */
    public static uploadLigandApiV1DrugDiscoveryUploadLigandPost(
        formData: Body_upload_ligand_api_v1_drug_discovery_upload_ligand_post,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/upload-ligand',
            formData: formData,
            mediaType: 'multipart/form-data',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Delete Ligand
     * @param experimentId
     * @param targetId
     * @param ligandId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static deleteLigandApiV1DrugDiscoveryDeleteLigandDelete(
        experimentId: string,
        targetId: string,
        ligandId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'DELETE',
            url: '/api/v1/drug_discovery/delete-ligand',
            query: {
                'experiment_id': experimentId,
                'target_id': targetId,
                'ligand_id': ligandId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get Ligands List
     * @param experimentId
     * @param targetId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static getLigandsListApiV1DrugDiscoveryGetLigandsListGet(
        experimentId: string,
        targetId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/get-ligands-list',
            query: {
                'experiment_id': experimentId,
                'target_id': targetId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get Target Binding Pocket
     * @param experimentId
     * @param targetId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static getTargetBindingPocketApiV1DrugDiscoveryGetTargetBindingPocketGet(
        experimentId: string,
        targetId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/get-target-binding-pocket',
            query: {
                'experiment_id': experimentId,
                'target_id': targetId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Predict Binding Pocket
     * @param requestBody
     * @returns any Successful Response
     * @throws ApiError
     */
    public static predictBindingPocketApiV1DrugDiscoveryPredictBindingPocketPost(
        requestBody: PredictBindingPocketRequest,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/predict-binding-pocket',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Predict Msa
     * @param requestBody
     * @returns any Successful Response
     * @throws ApiError
     */
    public static predictMsaApiV1DrugDiscoveryPredictMsaPost(
        requestBody: PredictMsaRequest,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/predict-msa',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Check Folding Exist
     * @param requestBody
     * @returns any Successful Response
     * @throws ApiError
     */
    public static checkFoldingExistApiV1DrugDiscoveryGetFoldedStructureGet(
        requestBody: GetFoldingRequest,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/get-folded-structure',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Predict Folding
     * @param requestBody
     * @returns any Successful Response
     * @throws ApiError
     */
    public static predictFoldingApiV1DrugDiscoveryPredictFoldingPost(
        requestBody: PredictFoldingRequest,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/predict-folding',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Perform Docking
     * @param requestBody
     * @returns any Successful Response
     * @throws ApiError
     */
    public static performDockingApiV1DrugDiscoveryPredictDockingPost(
        requestBody: DockingRequest,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/predict-docking',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
}
