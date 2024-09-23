/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { GetJobMetadataResponse } from '../models/GetJobMetadataResponse';
import type { UpdateJobRequest } from '../models/UpdateJobRequest';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class JobsACommonControllerForJobsManagementService {
    /**
     * Get all jobs metadata by experiment
     * @param experimentId
     * @returns GetJobMetadataResponse Successful Response
     * @throws ApiError
     */
    public static jobsMetadataApiV1JobsMetadataGet(
        experimentId: string,
    ): CancelablePromise<Array<GetJobMetadataResponse>> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/jobs/metadata',
            query: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get job metadata by experiment
     * @param jobId
     * @returns GetJobMetadataResponse Successful Response
     * @throws ApiError
     */
    public static jobMetadataApiV1JobsJobIdMetadataGet(
        jobId: string,
    ): CancelablePromise<GetJobMetadataResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/jobs/{job_id}/metadata',
            path: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Delete job
     * @param jobId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static deleteJobApiV1JobsJodIdDelete(
        jobId: string,
    ): CancelablePromise<any> {
      debugger;
        return __request(OpenAPI, {
            method: 'DELETE',
            url: '/api/v1/jobs/{jod_id}',
            query: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Update job
     * @param jobId
     * @param requestBody
     * @returns any Successful Response
     * @throws ApiError
     */
    public static updateApiV1JobsJobIdPatch(
        jobId: string,
        requestBody: UpdateJobRequest,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'PATCH',
            url: '/api/v1/jobs/{job_id}',
            path: {
                'job_id': jobId,
            },
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
}
