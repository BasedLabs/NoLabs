/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { application__folding__api_models__GetJobStatusResponse } from '../models/application__folding__api_models__GetJobStatusResponse';
import type { application__folding__api_models__JobResponse } from '../models/application__folding__api_models__JobResponse';
import type { application__folding__api_models__SetupJobRequest } from '../models/application__folding__api_models__SetupJobRequest';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class FoldingService {
    /**
     * Start folding job
     * @param jobId
     * @returns application__folding__api_models__JobResponse Successful Response
     * @throws ApiError
     */
    public static startJobApiV1FoldingJobsRunJobIdPost(
        jobId: string,
    ): CancelablePromise<application__folding__api_models__JobResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/folding/jobs/run/{job_id}',
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
     * @returns application__folding__api_models__JobResponse Successful Response
     * @throws ApiError
     */
    public static getJobApiV1FoldingJobsJobIdGet(
        jobId: string,
    ): CancelablePromise<application__folding__api_models__JobResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/folding/jobs/{job_id}',
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
     * @returns application__folding__api_models__GetJobStatusResponse Successful Response
     * @throws ApiError
     */
    public static getJobStatusApiV1FoldingJobsJobIdStatusGet(
        jobId: string,
    ): CancelablePromise<application__folding__api_models__GetJobStatusResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/folding/jobs/{job_id}/status',
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
     * @returns application__folding__api_models__JobResponse Successful Response
     * @throws ApiError
     */
    public static setupJobApiV1FoldingJobsPost(
        requestBody: application__folding__api_models__SetupJobRequest,
    ): CancelablePromise<application__folding__api_models__JobResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/folding/jobs',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
}
