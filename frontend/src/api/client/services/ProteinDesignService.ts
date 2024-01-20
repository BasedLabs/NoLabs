/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { Body_inference_api_v1_protein_design_inference_post } from '../models/Body_inference_api_v1_protein_design_inference_post';
import type { ChangeExperimentNameRequest } from '../models/ChangeExperimentNameRequest';
import type { GenerateUuidResponse } from '../models/GenerateUuidResponse';
import type { nolabs__api_models__protein_design__ExperimentMetadataResponse } from '../models/nolabs__api_models__protein_design__ExperimentMetadataResponse';
import type { nolabs__api_models__protein_design__GetExperimentResponse } from '../models/nolabs__api_models__protein_design__GetExperimentResponse';
import type { RunProteinDesignResponse } from '../models/RunProteinDesignResponse';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class ProteinDesignService {
    /**
     * Inference
     * @param formData
     * @returns RunProteinDesignResponse Successful Response
     * @throws ApiError
     */
    public static inferenceApiV1ProteinDesignInferencePost(
        formData: Body_inference_api_v1_protein_design_inference_post,
    ): CancelablePromise<RunProteinDesignResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/protein-design/inference',
            formData: formData,
            mediaType: 'multipart/form-data',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Experiments
     * @returns nolabs__api_models__protein_design__ExperimentMetadataResponse Successful Response
     * @throws ApiError
     */
    public static experimentsApiV1ProteinDesignExperimentsGet(): CancelablePromise<Array<nolabs__api_models__protein_design__ExperimentMetadataResponse>> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/protein-design/experiments',
        });
    }
    /**
     * Get Experiment
     * @param experimentId
     * @returns nolabs__api_models__protein_design__GetExperimentResponse Successful Response
     * @throws ApiError
     */
    public static getExperimentApiV1ProteinDesignGetExperimentGet(
        experimentId: string,
    ): CancelablePromise<nolabs__api_models__protein_design__GetExperimentResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/protein-design/get-experiment',
            query: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Delete Experiment
     * @param experimentId
     * @returns nolabs__api_models__protein_design__GetExperimentResponse Successful Response
     * @throws ApiError
     */
    public static deleteExperimentApiV1ProteinDesignDeleteExperimentDelete(
        experimentId: string,
    ): CancelablePromise<nolabs__api_models__protein_design__GetExperimentResponse> {
        return __request(OpenAPI, {
            method: 'DELETE',
            url: '/api/v1/protein-design/delete-experiment',
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
     * @returns nolabs__api_models__protein_design__GetExperimentResponse Successful Response
     * @throws ApiError
     */
    public static changeExperimentNameApiV1ProteinDesignChangeExperimentNamePost(
        requestBody: ChangeExperimentNameRequest,
    ): CancelablePromise<nolabs__api_models__protein_design__GetExperimentResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/protein-design/change-experiment-name',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Generate Uuid
     * @returns GenerateUuidResponse Successful Response
     * @throws ApiError
     */
    public static generateUuidApiV1ProteinDesignGenerateIdGet(): CancelablePromise<GenerateUuidResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/protein-design/generate_id',
        });
    }
}
