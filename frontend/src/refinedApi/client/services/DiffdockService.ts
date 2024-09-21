/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { application__diffdock__api_models__GetJobStatusResponse } from '../models/application__diffdock__api_models__GetJobStatusResponse';
import type { application__diffdock__api_models__JobResponse } from '../models/application__diffdock__api_models__JobResponse';
import type { application__diffdock__api_models__SetupJobRequest } from '../models/application__diffdock__api_models__SetupJobRequest';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class DiffdockService {
    /**
     * Start diffdock job
     * @param jobId
     * @returns application__diffdock__api_models__JobResponse Successful Response
     * @throws ApiError
     */
    public static startJobApiV1DiffdockJobsRunJobIdPost(
        jobId: string,
    ): CancelablePromise<application__diffdock__api_models__JobResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/diffdock/jobs/run/{job_id}',
            path: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get job
     * @param jobId
     * @returns application__diffdock__api_models__JobResponse Successful Response
     * @throws ApiError
     */
    public static getJobApiV1DiffdockJobsJobIdGet(
        jobId: string,
    ): CancelablePromise<application__diffdock__api_models__JobResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/diffdock/jobs/{job_id}',
            path: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get job execution status
     * @param jobId
     * @returns application__diffdock__api_models__GetJobStatusResponse Successful Response
     * @throws ApiError
     */
    public static getJobStatusApiV1DiffdockJobsJobIdStatusGet(
        jobId: string,
    ): CancelablePromise<application__diffdock__api_models__GetJobStatusResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/diffdock/jobs/{job_id}/status',
            path: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Setup job
     * @param requestBody
     * @returns application__diffdock__api_models__JobResponse Successful Response
     * @throws ApiError
     */
    public static setupJobApiV1DiffdockJobsPost(
        requestBody: application__diffdock__api_models__SetupJobRequest,
    ): CancelablePromise<application__diffdock__api_models__JobResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/diffdock/jobs',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
}
