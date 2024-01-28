/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { Body_upload_ligand_api_v1_drug_discovery_upload_ligand_post } from '../models/Body_upload_ligand_api_v1_drug_discovery_upload_ligand_post';
import type { Body_upload_target_api_v1_drug_discovery_upload_target_post } from '../models/Body_upload_target_api_v1_drug_discovery_upload_target_post';
import type { ChangeExperimentNameRequest } from '../models/ChangeExperimentNameRequest';
import type { DeleteLigandResponse } from '../models/DeleteLigandResponse';
import type { DeleteTargetResponse } from '../models/DeleteTargetResponse';
import type { DockingRequest } from '../models/DockingRequest';
import type { DockingResponse } from '../models/DockingResponse';
import type { ExperimentMetadataResponse } from '../models/ExperimentMetadataResponse';
import type { GetFoldingRequest } from '../models/GetFoldingRequest';
import type { GetFoldingResponse } from '../models/GetFoldingResponse';
import type { GetLigandDataResponse } from '../models/GetLigandDataResponse';
import type { GetTargetBindingPocketResponse } from '../models/GetTargetBindingPocketResponse';
import type { GetTargetDataResponse } from '../models/GetTargetDataResponse';
import type { LigandMetaData } from '../models/LigandMetaData';
import type { PredictBindingPocketResponse } from '../models/PredictBindingPocketResponse';
import type { PredictFoldingResponse } from '../models/PredictFoldingResponse';
import type { PredictMsaResponse } from '../models/PredictMsaResponse';
import type { TargetMetaData } from '../models/TargetMetaData';
import type { UploadLigandResponse } from '../models/UploadLigandResponse';
import type { UploadTargetResponse } from '../models/UploadTargetResponse';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class DrugDiscoveryService {
    /**
     * Experiments
     * @returns ExperimentMetadataResponse Successful Response
     * @throws ApiError
     */
    public static experimentsApiV1DrugDiscoveryExperimentsGet(): CancelablePromise<Array<ExperimentMetadataResponse>> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/experiments',
        });
    }
    /**
     * Add Experiment
     * @returns ExperimentMetadataResponse Successful Response
     * @throws ApiError
     */
    public static addExperimentApiV1DrugDiscoveryAddExperimentPost(): CancelablePromise<ExperimentMetadataResponse> {
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
     * @returns ExperimentMetadataResponse Successful Response
     * @throws ApiError
     */
    public static changeExperimentNameApiV1DrugDiscoveryChangeExperimentNamePost(
        requestBody: ChangeExperimentNameRequest,
    ): CancelablePromise<ExperimentMetadataResponse> {
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
     * @returns UploadTargetResponse Successful Response
     * @throws ApiError
     */
    public static uploadTargetApiV1DrugDiscoveryUploadTargetPost(
        formData: Body_upload_target_api_v1_drug_discovery_upload_target_post,
    ): CancelablePromise<UploadTargetResponse> {
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
     * @returns DeleteTargetResponse Successful Response
     * @throws ApiError
     */
    public static deleteTargetApiV1DrugDiscoveryDeleteTargetDelete(
        experimentId: string,
        targetId: string,
    ): CancelablePromise<DeleteTargetResponse> {
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
     * @returns TargetMetaData Successful Response
     * @throws ApiError
     */
    public static getTargetsListApiV1DrugDiscoveryGetTargetsListGet(
        experimentId: string,
    ): CancelablePromise<Array<TargetMetaData>> {
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
     * Get Targets List
     * @param experimentId
     * @param targetId
     * @returns GetTargetDataResponse Successful Response
     * @throws ApiError
     */
    public static getTargetsListApiV1DrugDiscoveryGetTargetDataGet(
        experimentId: string,
        targetId: string,
    ): CancelablePromise<GetTargetDataResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/get-target-data',
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
     * Upload Ligand
     * @param formData
     * @returns UploadLigandResponse Successful Response
     * @throws ApiError
     */
    public static uploadLigandApiV1DrugDiscoveryUploadLigandPost(
        formData: Body_upload_ligand_api_v1_drug_discovery_upload_ligand_post,
    ): CancelablePromise<UploadLigandResponse> {
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
     * @returns DeleteLigandResponse Successful Response
     * @throws ApiError
     */
    public static deleteLigandApiV1DrugDiscoveryDeleteLigandDelete(
        experimentId: string,
        targetId: string,
        ligandId: string,
    ): CancelablePromise<DeleteLigandResponse> {
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
     * Get Ligand Data
     * @param experimentId
     * @param targetId
     * @param ligandId
     * @returns GetLigandDataResponse Successful Response
     * @throws ApiError
     */
    public static getLigandDataApiV1DrugDiscoveryGetLigandDataGet(
        experimentId: string,
        targetId: string,
        ligandId: string,
    ): CancelablePromise<GetLigandDataResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/get-ligand-data',
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
     * @returns LigandMetaData Successful Response
     * @throws ApiError
     */
    public static getLigandsListApiV1DrugDiscoveryGetLigandsListGet(
        experimentId: string,
        targetId: string,
    ): CancelablePromise<Array<LigandMetaData>> {
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
     * @returns GetTargetBindingPocketResponse Successful Response
     * @throws ApiError
     */
    public static getTargetBindingPocketApiV1DrugDiscoveryGetTargetBindingPocketGet(
        experimentId: string,
        targetId: string,
    ): CancelablePromise<GetTargetBindingPocketResponse> {
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
     * @param experimentId
     * @param targetId
     * @returns PredictBindingPocketResponse Successful Response
     * @throws ApiError
     */
    public static predictBindingPocketApiV1DrugDiscoveryPredictBindingPocketPost(
        experimentId: string,
        targetId: string,
    ): CancelablePromise<PredictBindingPocketResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/predict-binding-pocket',
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
     * Predict Msa
     * @param experimentId
     * @param targetId
     * @returns PredictMsaResponse Successful Response
     * @throws ApiError
     */
    public static predictMsaApiV1DrugDiscoveryPredictMsaPost(
        experimentId: string,
        targetId: string,
    ): CancelablePromise<PredictMsaResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/predict-msa',
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
     * Get Folded Structure
     * @param requestBody
     * @returns GetFoldingResponse Successful Response
     * @throws ApiError
     */
    public static getFoldedStructureApiV1DrugDiscoveryGetFoldedStructureGet(
        requestBody: GetFoldingRequest,
    ): CancelablePromise<GetFoldingResponse> {
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
     * @param experimentId
     * @param targetId
     * @returns PredictFoldingResponse Successful Response
     * @throws ApiError
     */
    public static predictFoldingApiV1DrugDiscoveryPredictFoldingPost(
        experimentId: string,
        targetId: string,
    ): CancelablePromise<PredictFoldingResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/predict-folding',
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
     * Perform Docking
     * @param requestBody
     * @returns DockingResponse Successful Response
     * @throws ApiError
     */
    public static performDockingApiV1DrugDiscoveryPredictDockingPost(
        requestBody: DockingRequest,
    ): CancelablePromise<DockingResponse> {
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
