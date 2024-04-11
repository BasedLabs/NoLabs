/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { Body_inference_api_v1_folding_inference_post } from '../models/Body_inference_api_v1_folding_inference_post';
import type { ChangeExperimentNameRequest } from '../models/ChangeExperimentNameRequest';
import type { ExperimentMetadataResponse } from '../models/ExperimentMetadataResponse';
import type { nolabs__api_models__amino_acid__folding__GetExperimentResponse } from '../models/nolabs__api_models__amino_acid__folding__GetExperimentResponse';
import type { RunFoldingResponse } from '../models/RunFoldingResponse';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class FoldingService {
    /**
     * Inference
     * @param formData
     * @returns RunFoldingResponse Successful Response
     * @throws ApiError
     */
    public static inferenceApiV1FoldingInferencePost(
        formData: Body_inference_api_v1_folding_inference_post,
    ): CancelablePromise<RunFoldingResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/folding/inference',
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
    public static experimentsApiV1FoldingExperimentsGet(): CancelablePromise<Array<ExperimentMetadataResponse>> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/folding/experiments',
        });
    }
    /**
     * Get Experiment
     * @param experimentId
     * @returns nolabs__api_models__amino_acid__folding__GetExperimentResponse Successful Response
     * @throws ApiError
     */
    public static getExperimentApiV1FoldingGetExperimentGet(
        experimentId: string,
    ): CancelablePromise<nolabs__api_models__amino_acid__folding__GetExperimentResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/folding/get-experiment',
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
     * @returns any Successful Response
     * @throws ApiError
     */
    public static deleteExperimentApiV1FoldingDeleteExperimentDelete(
        experimentId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'DELETE',
            url: '/api/v1/folding/delete-experiment',
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
    public static changeExperimentNameApiV1FoldingChangeExperimentNamePost(
        requestBody: ChangeExperimentNameRequest,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/folding/change-experiment-name',
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
    public static createExperimentApiV1FoldingCreateExperimentGet(): CancelablePromise<ExperimentMetadataResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/folding/create-experiment',
        });
    }
}
