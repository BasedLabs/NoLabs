/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { Body_inference_api_v1_protein_design_inference_post } from '../models/Body_inference_api_v1_protein_design_inference_post';
import type { ChangeExperimentNameRequest } from '../models/ChangeExperimentNameRequest';
import type { ExperimentMetadataResponse } from '../models/ExperimentMetadataResponse';
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
     * @returns ExperimentMetadataResponse Successful Response
     * @throws ApiError
     */
    public static experimentsApiV1ProteinDesignExperimentsMetadataGet(): CancelablePromise<Array<ExperimentMetadataResponse>> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/protein-design/experiments-metadata',
        });
    }
    /**
     * Get Experiment
     * @param id
     * @returns nolabs__api_models__protein_design__GetExperimentResponse Successful Response
     * @throws ApiError
     */
    public static getExperimentApiV1ProteinDesignExperimentGet(
        id: string,
    ): CancelablePromise<nolabs__api_models__protein_design__GetExperimentResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/protein-design/experiment',
            query: {
                'id': id,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Delete Experiment
     * @param id
     * @returns nolabs__api_models__protein_design__GetExperimentResponse Successful Response
     * @throws ApiError
     */
    public static deleteExperimentApiV1ProteinDesignExperimentDelete(
        id: string,
    ): CancelablePromise<nolabs__api_models__protein_design__GetExperimentResponse> {
        return __request(OpenAPI, {
            method: 'DELETE',
            url: '/api/v1/protein-design/experiment',
            query: {
                'id': id,
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
    public static changeExperimentNameApiV1ProteinDesignChangeExperimentNamePost(
        requestBody: ChangeExperimentNameRequest,
    ): CancelablePromise<any> {
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
     * Create Experiment
     * @returns ExperimentMetadataResponse Successful Response
     * @throws ApiError
     */
    public static createExperimentApiV1ProteinDesignCreateExperimentGet(): CancelablePromise<ExperimentMetadataResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/protein-design/create-experiment',
        });
    }
}
