/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { Body_inference_api_v1_localisation_inference_post } from '../models/Body_inference_api_v1_localisation_inference_post';
import type { ChangeExperimentNameRequest } from '../models/ChangeExperimentNameRequest';
import type { ExperimentMetadataResponse } from '../models/ExperimentMetadataResponse';
import type { nolabs__api_models__amino_acid__localisation__GetExperimentResponse } from '../models/nolabs__api_models__amino_acid__localisation__GetExperimentResponse';
import type { RunLocalisationResponse } from '../models/RunLocalisationResponse';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class LocalisationService {
    /**
     * Inference
     * @param formData
     * @returns RunLocalisationResponse Successful Response
     * @throws ApiError
     */
    public static inferenceApiV1LocalisationInferencePost(
        formData: Body_inference_api_v1_localisation_inference_post,
    ): CancelablePromise<RunLocalisationResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/localisation/inference',
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
    public static experimentsApiV1LocalisationExperimentsGet(): CancelablePromise<Array<ExperimentMetadataResponse>> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/localisation/experiments',
        });
    }
    /**
     * Get Experiment
     * @param experimentId
     * @returns nolabs__api_models__amino_acid__localisation__GetExperimentResponse Successful Response
     * @throws ApiError
     */
    public static getExperimentApiV1LocalisationGetExperimentGet(
        experimentId: string,
    ): CancelablePromise<nolabs__api_models__amino_acid__localisation__GetExperimentResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/localisation/get-experiment',
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
    public static deleteExperimentApiV1LocalisationDeleteExperimentDelete(
        experimentId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'DELETE',
            url: '/api/v1/localisation/delete-experiment',
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
    public static changeExperimentNameApiV1LocalisationChangeExperimentNamePost(
        requestBody: ChangeExperimentNameRequest,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/localisation/change-experiment-name',
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
    public static createExperimentApiV1LocalisationCreateExperimentGet(): CancelablePromise<ExperimentMetadataResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/localisation/create-experiment',
        });
    }
}
