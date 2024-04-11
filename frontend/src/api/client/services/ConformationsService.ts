/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { Body_inference_api_v1_conformations_inference_post } from '../models/Body_inference_api_v1_conformations_inference_post';
import type { ChangeExperimentNameRequest } from '../models/ChangeExperimentNameRequest';
import type { ExperimentMetadataResponse } from '../models/ExperimentMetadataResponse';
import type { nolabs__api_models__conformations__GetExperimentResponse } from '../models/nolabs__api_models__conformations__GetExperimentResponse';
import type { RunSimulationsResponse } from '../models/RunSimulationsResponse';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class ConformationsService {
    /**
     * Inference
     * @param formData
     * @returns RunSimulationsResponse Successful Response
     * @throws ApiError
     */
    public static inferenceApiV1ConformationsInferencePost(
        formData: Body_inference_api_v1_conformations_inference_post,
    ): CancelablePromise<RunSimulationsResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/conformations/inference',
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
    public static experimentsApiV1ConformationsExperimentsGet(): CancelablePromise<Array<ExperimentMetadataResponse>> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/conformations/experiments',
        });
    }
    /**
     * Get Experiment
     * @param experimentId
     * @returns nolabs__api_models__conformations__GetExperimentResponse Successful Response
     * @throws ApiError
     */
    public static getExperimentApiV1ConformationsGetExperimentGet(
        experimentId: string,
    ): CancelablePromise<nolabs__api_models__conformations__GetExperimentResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/conformations/get-experiment',
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
    public static deleteExperimentApiV1ConformationsDeleteExperimentDelete(
        experimentId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'DELETE',
            url: '/api/v1/conformations/delete-experiment',
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
    public static changeExperimentNameApiV1ConformationsChangeExperimentNamePost(
        requestBody: ChangeExperimentNameRequest,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/conformations/change-experiment-name',
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
    public static createExperimentApiV1ConformationsCreateExperimentGet(): CancelablePromise<ExperimentMetadataResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/conformations/create-experiment',
        });
    }
}
