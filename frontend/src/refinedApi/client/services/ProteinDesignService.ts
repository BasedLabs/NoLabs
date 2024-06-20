/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { nolabs__application__use_cases__protein_design__api_models__JobResponse } from '../models/nolabs__application__use_cases__protein_design__api_models__JobResponse';
import type { nolabs__application__use_cases__protein_design__api_models__SetupJobRequest } from '../models/nolabs__application__use_cases__protein_design__api_models__SetupJobRequest';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class ProteinDesignService {
    /**
     * Run Job
     * @param jobId
     * @returns nolabs__application__use_cases__protein_design__api_models__JobResponse Successful Response
     * @throws ApiError
     */
    public static runJobApiV1ProteinDesignJobsRunJobIdPost(
        jobId: string,
    ): CancelablePromise<nolabs__application__use_cases__protein_design__api_models__JobResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/protein-design/jobs/run/{job_id}',
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
     * @returns nolabs__application__use_cases__protein_design__api_models__JobResponse Successful Response
     * @throws ApiError
     */
    public static getJobApiV1ProteinDesignJobsJobIdGet(
        jobId: string,
    ): CancelablePromise<nolabs__application__use_cases__protein_design__api_models__JobResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/protein-design/jobs/{job_id}',
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
     * @returns nolabs__application__use_cases__protein_design__api_models__JobResponse Successful Response
     * @throws ApiError
     */
    public static setupJobApiV1ProteinDesignJobsPost(
        requestBody: nolabs__application__use_cases__protein_design__api_models__SetupJobRequest,
    ): CancelablePromise<nolabs__application__use_cases__protein_design__api_models__JobResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/protein-design/jobs',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
}
