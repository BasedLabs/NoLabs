/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { Body_upload_ligand_api_v1_drug_discovery_upload_ligand_post } from '../models/Body_upload_ligand_api_v1_drug_discovery_upload_ligand_post';
import type { Body_upload_target_api_v1_drug_discovery_upload_target_post } from '../models/Body_upload_target_api_v1_drug_discovery_upload_target_post';
import type { CheckFoldingDataAvailableResponse } from '../models/CheckFoldingDataAvailableResponse';
import type { CheckJobIsRunningResponse } from '../models/CheckJobIsRunningResponse';
import type { CheckMsaDataAvailableResponse } from '../models/CheckMsaDataAvailableResponse';
import type { CheckPocketDataAvailableResponse } from '../models/CheckPocketDataAvailableResponse';
import type { CheckResultDataAvailableResponse } from '../models/CheckResultDataAvailableResponse';
import type { CheckServiceHealthyResponse } from '../models/CheckServiceHealthyResponse';
import type { DeleteDockingJobResponse } from '../models/DeleteDockingJobResponse';
import type { DeleteLigandResponse } from '../models/DeleteLigandResponse';
import type { DeleteTargetResponse } from '../models/DeleteTargetResponse';
import type { ExperimentMetadataResponse } from '../models/ExperimentMetadataResponse';
import type { GetAllJobsListResponse } from '../models/GetAllJobsListResponse';
import type { GetAllResultsListResponse } from '../models/GetAllResultsListResponse';
import type { GetDiffDockDockingResultDataResponse } from '../models/GetDiffDockDockingResultDataResponse';
import type { GetDiffDockLigandSdfResponse } from '../models/GetDiffDockLigandSdfResponse';
import type { GetFoldingResponse } from '../models/GetFoldingResponse';
import type { GetJobBindingPocketDataResponse } from '../models/GetJobBindingPocketDataResponse';
import type { GetJobsListForTargetLigandResponse } from '../models/GetJobsListForTargetLigandResponse';
import type { GetLigandDataResponse } from '../models/GetLigandDataResponse';
import type { GetLigandMetaDataResponse } from '../models/GetLigandMetaDataResponse';
import type { GetResultsListForTargetLigandResponse } from '../models/GetResultsListForTargetLigandResponse';
import type { GetTargetBindingPocketResponse } from '../models/GetTargetBindingPocketResponse';
import type { GetTargetDataResponse } from '../models/GetTargetDataResponse';
import type { GetTargetMetaDataResponse } from '../models/GetTargetMetaDataResponse';
import type { GetUmolDockingResultDataResponse } from '../models/GetUmolDockingResultDataResponse';
import type { LigandMetaData } from '../models/LigandMetaData';
import type { PredictBindingPocketResponse } from '../models/PredictBindingPocketResponse';
import type { PredictFoldingResponse } from '../models/PredictFoldingResponse';
import type { PredictMsaResponse } from '../models/PredictMsaResponse';
import type { RegisterDockingJobResponse } from '../models/RegisterDockingJobResponse';
import type { RunDiffDockDockingJobResponse } from '../models/RunDiffDockDockingJobResponse';
import type { RunUmolDockingJobResponse } from '../models/RunUmolDockingJobResponse';
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
     * Create Experiment
     * @returns ExperimentMetadataResponse Successful Response
     * @throws ApiError
     */
    public static createExperimentApiV1DrugDiscoveryCreateExperimentGet(): CancelablePromise<ExperimentMetadataResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/create-experiment',
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
     * Get Experiment Metadata
     * @param experimentId
     * @returns ExperimentMetadataResponse Successful Response
     * @throws ApiError
     */
    public static getExperimentMetadataApiV1DrugDiscoveryGetExperimentMetadataGet(
        experimentId: string,
    ): CancelablePromise<ExperimentMetadataResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/get-experiment-metadata',
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
     * @param experimentId
     * @param experimentName
     * @returns ExperimentMetadataResponse Successful Response
     * @throws ApiError
     */
    public static changeExperimentNameApiV1DrugDiscoveryChangeExperimentNamePost(
        experimentId: string,
        experimentName: string,
    ): CancelablePromise<ExperimentMetadataResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/change-experiment-name',
            query: {
                'experiment_id': experimentId,
                'experiment_name': experimentName,
            },
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
     * Get Target Meta Data
     * @param experimentId
     * @param targetId
     * @returns GetTargetMetaDataResponse Successful Response
     * @throws ApiError
     */
    public static getTargetMetaDataApiV1DrugDiscoveryGetTargetMetaDataGet(
        experimentId: string,
        targetId: string,
    ): CancelablePromise<GetTargetMetaDataResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/get-target-meta-data',
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
     * Update Target Name
     * @param experimentId
     * @param targetId
     * @param targetName
     * @returns GetTargetMetaDataResponse Successful Response
     * @throws ApiError
     */
    public static updateTargetNameApiV1DrugDiscoveryUpdateTargetNamePost(
        experimentId: string,
        targetId: string,
        targetName: string,
    ): CancelablePromise<GetTargetMetaDataResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/update-target-name',
            query: {
                'experiment_id': experimentId,
                'target_id': targetId,
                'target_name': targetName,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get Target Data
     * @param experimentId
     * @param targetId
     * @returns GetTargetDataResponse Successful Response
     * @throws ApiError
     */
    public static getTargetDataApiV1DrugDiscoveryGetTargetDataGet(
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
     * @returns GetLigandMetaDataResponse Successful Response
     * @throws ApiError
     */
    public static getLigandDataApiV1DrugDiscoveryGetLigandMetaDataGet(
        experimentId: string,
        targetId: string,
        ligandId: string,
    ): CancelablePromise<GetLigandMetaDataResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/get-ligand-meta-data',
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
     * Set Target Binding Pocket
     * @param experimentId
     * @param targetId
     * @param requestBody
     * @returns any Successful Response
     * @throws ApiError
     */
    public static setTargetBindingPocketApiV1DrugDiscoverySetTargetBindingPocketPost(
        experimentId: string,
        targetId: string,
        requestBody: Array<number>,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/set-target-binding-pocket',
            query: {
                'experiment_id': experimentId,
                'target_id': targetId,
            },
            body: requestBody,
            mediaType: 'application/json',
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
     * @param experimentId
     * @param targetId
     * @returns GetFoldingResponse Successful Response
     * @throws ApiError
     */
    public static getFoldedStructureApiV1DrugDiscoveryGetFoldedStructureGet(
        experimentId: string,
        targetId: string,
    ): CancelablePromise<GetFoldingResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/get-folded-structure',
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
     * Predict Folding
     * @param experimentId
     * @param targetId
     * @returns PredictFoldingResponse Successful Response
     * @throws ApiError
     */
    public static predictFoldingApiV1DrugDiscoveryPredictEsmfoldLightPost(
        experimentId: string,
        targetId: string,
    ): CancelablePromise<PredictFoldingResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/predict-esmfold-light',
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
     * Predict Folding
     * @param experimentId
     * @param targetId
     * @returns PredictFoldingResponse Successful Response
     * @throws ApiError
     */
    public static predictFoldingApiV1DrugDiscoveryPredictEsmfoldPost(
        experimentId: string,
        targetId: string,
    ): CancelablePromise<PredictFoldingResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/predict-esmfold',
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
     * Register Docking Job
     * @param experimentId
     * @param targetId
     * @param ligandId
     * @param foldingMethod
     * @returns RegisterDockingJobResponse Successful Response
     * @throws ApiError
     */
    public static registerDockingJobApiV1DrugDiscoveryRegisterDockingJobPost(
        experimentId: string,
        targetId: string,
        ligandId: string,
        foldingMethod: string,
    ): CancelablePromise<RegisterDockingJobResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/register-docking-job',
            query: {
                'experiment_id': experimentId,
                'target_id': targetId,
                'ligand_id': ligandId,
                'folding_method': foldingMethod,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Delete Docking Job
     * @param experimentId
     * @param targetId
     * @param ligandId
     * @param jobId
     * @returns DeleteDockingJobResponse Successful Response
     * @throws ApiError
     */
    public static deleteDockingJobApiV1DrugDiscoveryDeleteDockingJobDelete(
        experimentId: string,
        targetId: string,
        ligandId: string,
        jobId: string,
    ): CancelablePromise<DeleteDockingJobResponse> {
        return __request(OpenAPI, {
            method: 'DELETE',
            url: '/api/v1/drug_discovery/delete-docking-job',
            query: {
                'experiment_id': experimentId,
                'target_id': targetId,
                'ligand_id': ligandId,
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Check Msa Data Available
     * @param experimentId
     * @param targetId
     * @returns CheckMsaDataAvailableResponse Successful Response
     * @throws ApiError
     */
    public static checkMsaDataAvailableApiV1DrugDiscoveryCheckMsaDataAvailableGet(
        experimentId: string,
        targetId: string,
    ): CancelablePromise<CheckMsaDataAvailableResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/check-msa-data-available',
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
     * Check Msa Service Health
     * @returns CheckServiceHealthyResponse Successful Response
     * @throws ApiError
     */
    public static checkMsaServiceHealthApiV1DrugDiscoveryCheckMsaServiceHealthGet(): CancelablePromise<CheckServiceHealthyResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/check-msa-service-health',
        });
    }
    /**
     * Check Msa Job Running
     * @param jobId
     * @returns CheckJobIsRunningResponse Successful Response
     * @throws ApiError
     */
    public static checkMsaJobRunningApiV1DrugDiscoveryCheckMsaJobRunningGet(
        jobId: string,
    ): CancelablePromise<CheckJobIsRunningResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/check-msa-job-running',
            query: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Check Pocket Data Available
     * @param experimentId
     * @param targetId
     * @param ligandId
     * @param jobId
     * @returns CheckPocketDataAvailableResponse Successful Response
     * @throws ApiError
     */
    public static checkPocketDataAvailableApiV1DrugDiscoveryCheckPocketDataAvailableGet(
        experimentId: string,
        targetId: string,
        ligandId: string,
        jobId: string,
    ): CancelablePromise<CheckPocketDataAvailableResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/check-pocket-data-available',
            query: {
                'experiment_id': experimentId,
                'target_id': targetId,
                'ligand_id': ligandId,
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get Job Pocket Data
     * @param experimentId
     * @param targetId
     * @param ligandId
     * @param jobId
     * @returns GetJobBindingPocketDataResponse Successful Response
     * @throws ApiError
     */
    public static getJobPocketDataApiV1DrugDiscoveryGetJobPocketDataGet(
        experimentId: string,
        targetId: string,
        ligandId: string,
        jobId: string,
    ): CancelablePromise<GetJobBindingPocketDataResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/get-job-pocket-data',
            query: {
                'experiment_id': experimentId,
                'target_id': targetId,
                'ligand_id': ligandId,
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Check P2Rank Service Health
     * @returns CheckServiceHealthyResponse Successful Response
     * @throws ApiError
     */
    public static checkP2RankServiceHealthApiV1DrugDiscoveryCheckP2RankServiceHealthGet(): CancelablePromise<CheckServiceHealthyResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/check-p2rank-service-health',
        });
    }
    /**
     * Check P2Rank Job Running
     * @param jobId
     * @returns CheckJobIsRunningResponse Successful Response
     * @throws ApiError
     */
    public static checkP2RankJobRunningApiV1DrugDiscoveryCheckP2RankJobRunningGet(
        jobId: string,
    ): CancelablePromise<CheckJobIsRunningResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/check-p2rank-job-running',
            query: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Check Folding Data Available
     * @param experimentId
     * @param targetId
     * @param foldingMethod
     * @returns CheckFoldingDataAvailableResponse Successful Response
     * @throws ApiError
     */
    public static checkFoldingDataAvailableApiV1DrugDiscoveryCheckFoldingDataAvailableGet(
        experimentId: string,
        targetId: string,
        foldingMethod: string,
    ): CancelablePromise<CheckFoldingDataAvailableResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/check-folding-data-available',
            query: {
                'experiment_id': experimentId,
                'target_id': targetId,
                'folding_method': foldingMethod,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Check Esmfold Service Health
     * @returns CheckServiceHealthyResponse Successful Response
     * @throws ApiError
     */
    public static checkEsmfoldServiceHealthApiV1DrugDiscoveryCheckEsmfoldServiceHealthGet(): CancelablePromise<CheckServiceHealthyResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/check-esmfold-service-health',
        });
    }
    /**
     * Check Esmfold Light Service Health
     * @returns CheckServiceHealthyResponse Successful Response
     * @throws ApiError
     */
    public static checkEsmfoldLightServiceHealthApiV1DrugDiscoveryCheckEsmfoldLightServiceHealthGet(): CancelablePromise<CheckServiceHealthyResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/check-esmfold-light-service-health',
        });
    }
    /**
     * Check Esmfold Job Running
     * @param jobId
     * @returns CheckJobIsRunningResponse Successful Response
     * @throws ApiError
     */
    public static checkEsmfoldJobRunningApiV1DrugDiscoveryCheckEsmfoldJobRunningGet(
        jobId: string,
    ): CancelablePromise<CheckJobIsRunningResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/check-esmfold-job-running',
            query: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Check Esmfold Light Job Running
     * @param jobId
     * @returns CheckJobIsRunningResponse Successful Response
     * @throws ApiError
     */
    public static checkEsmfoldLightJobRunningApiV1DrugDiscoveryCheckEsmfoldLightJobRunningGet(
        jobId: string,
    ): CancelablePromise<CheckJobIsRunningResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/check-esmfold-light-job-running',
            query: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Check Umol Service Health
     * @returns CheckServiceHealthyResponse Successful Response
     * @throws ApiError
     */
    public static checkUmolServiceHealthApiV1DrugDiscoveryCheckUmolServiceHealthGet(): CancelablePromise<CheckServiceHealthyResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/check-umol-service-health',
        });
    }
    /**
     * Check Umol Job Running
     * @param jobId
     * @returns CheckJobIsRunningResponse Successful Response
     * @throws ApiError
     */
    public static checkUmolJobRunningApiV1DrugDiscoveryCheckUmolJobRunningGet(
        jobId: string,
    ): CancelablePromise<CheckJobIsRunningResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/check-umol-job-running',
            query: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Check Diffdock Service Health
     * @returns CheckServiceHealthyResponse Successful Response
     * @throws ApiError
     */
    public static checkDiffdockServiceHealthApiV1DrugDiscoveryCheckDiffdockServiceHealthGet(): CancelablePromise<CheckServiceHealthyResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/check-diffdock-service-health',
        });
    }
    /**
     * Check Diffdock Job Running
     * @param jobId
     * @returns CheckJobIsRunningResponse Successful Response
     * @throws ApiError
     */
    public static checkDiffdockJobRunningApiV1DrugDiscoveryCheckDiffdockJobRunningGet(
        jobId: string,
    ): CancelablePromise<CheckJobIsRunningResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/check-diffdock-job-running',
            query: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Check Result Data Available
     * @param experimentId
     * @param targetId
     * @param ligandId
     * @param jobId
     * @returns CheckResultDataAvailableResponse Successful Response
     * @throws ApiError
     */
    public static checkResultDataAvailableApiV1DrugDiscoveryCheckResultDataAvailableGet(
        experimentId: string,
        targetId: string,
        ligandId: string,
        jobId: string,
    ): CancelablePromise<CheckResultDataAvailableResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/check-result-data-available',
            query: {
                'experiment_id': experimentId,
                'target_id': targetId,
                'ligand_id': ligandId,
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Run Umol Docking
     * @param experimentId
     * @param targetId
     * @param ligandId
     * @param jobId
     * @returns RunUmolDockingJobResponse Successful Response
     * @throws ApiError
     */
    public static runUmolDockingApiV1DrugDiscoveryRunUmolDockingJobPost(
        experimentId: string,
        targetId: string,
        ligandId: string,
        jobId: string,
    ): CancelablePromise<RunUmolDockingJobResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/run-umol-docking-job',
            query: {
                'experiment_id': experimentId,
                'target_id': targetId,
                'ligand_id': ligandId,
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Run Diffdock Docking
     * @param experimentId
     * @param targetId
     * @param ligandId
     * @param jobId
     * @returns RunDiffDockDockingJobResponse Successful Response
     * @throws ApiError
     */
    public static runDiffdockDockingApiV1DrugDiscoveryRunDiffdockDockingJobPost(
        experimentId: string,
        targetId: string,
        ligandId: string,
        jobId: string,
    ): CancelablePromise<RunDiffDockDockingJobResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/drug_discovery/run-diffdock-docking-job',
            query: {
                'experiment_id': experimentId,
                'target_id': targetId,
                'ligand_id': ligandId,
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get Results List For Target Ligand
     * @param experimentId
     * @param targetId
     * @param ligandId
     * @returns GetResultsListForTargetLigandResponse Successful Response
     * @throws ApiError
     */
    public static getResultsListForTargetLigandApiV1DrugDiscoveryGetResultsListForTargetLigandGet(
        experimentId: string,
        targetId: string,
        ligandId: string,
    ): CancelablePromise<GetResultsListForTargetLigandResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/get-results-list-for-target-ligand',
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
     * Get Jobs List For Target Ligand
     * @param experimentId
     * @param targetId
     * @param ligandId
     * @returns GetJobsListForTargetLigandResponse Successful Response
     * @throws ApiError
     */
    public static getJobsListForTargetLigandApiV1DrugDiscoveryGetJobsListForTargetLigandGet(
        experimentId: string,
        targetId: string,
        ligandId: string,
    ): CancelablePromise<GetJobsListForTargetLigandResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/get-jobs-list-for-target-ligand',
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
     * Get All Results List
     * @param experimentId
     * @returns GetAllResultsListResponse Successful Response
     * @throws ApiError
     */
    public static getAllResultsListApiV1DrugDiscoveryGetAllResultsListGet(
        experimentId: string,
    ): CancelablePromise<GetAllResultsListResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/get-all-results-list',
            query: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get All Jobs List
     * @param experimentId
     * @returns GetAllJobsListResponse Successful Response
     * @throws ApiError
     */
    public static getAllJobsListApiV1DrugDiscoveryGetAllJobsListGet(
        experimentId: string,
    ): CancelablePromise<GetAllJobsListResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/get-all-jobs-list',
            query: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get Umol Docking Result Data
     * @param experimentId
     * @param targetId
     * @param ligandId
     * @param jobId
     * @returns GetUmolDockingResultDataResponse Successful Response
     * @throws ApiError
     */
    public static getUmolDockingResultDataApiV1DrugDiscoveryGetUmolDockingResultDataGet(
        experimentId: string,
        targetId: string,
        ligandId: string,
        jobId: string,
    ): CancelablePromise<GetUmolDockingResultDataResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/get-umol-docking-result-data',
            query: {
                'experiment_id': experimentId,
                'target_id': targetId,
                'ligand_id': ligandId,
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get Diffdock Docking Result Data
     * @param experimentId
     * @param targetId
     * @param ligandId
     * @param jobId
     * @returns GetDiffDockDockingResultDataResponse Successful Response
     * @throws ApiError
     */
    public static getDiffdockDockingResultDataApiV1DrugDiscoveryGetDiffdockDockingResultDataGet(
        experimentId: string,
        targetId: string,
        ligandId: string,
        jobId: string,
    ): CancelablePromise<GetDiffDockDockingResultDataResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/get-diffdock-docking-result-data',
            query: {
                'experiment_id': experimentId,
                'target_id': targetId,
                'ligand_id': ligandId,
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get Diffdock Ligand Sdf
     * @param experimentId
     * @param targetId
     * @param ligandId
     * @param jobId
     * @param sdfFileName
     * @returns GetDiffDockLigandSdfResponse Successful Response
     * @throws ApiError
     */
    public static getDiffdockLigandSdfApiV1DrugDiscoveryGetDiffdockLigandSdfGet(
        experimentId: string,
        targetId: string,
        ligandId: string,
        jobId: string,
        sdfFileName: string,
    ): CancelablePromise<GetDiffDockLigandSdfResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/drug_discovery/get-diffdock-ligand_sdf',
            query: {
                'experiment_id': experimentId,
                'target_id': targetId,
                'ligand_id': ligandId,
                'job_id': jobId,
                'sdf_file_name': sdfFileName,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
}
