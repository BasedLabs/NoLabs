/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { EstimatesResponse } from '../models/EstimatesResponse';
import type { nolabs__application__adaptyv_bio__protein_binding_screening__api_models__JobResponse } from '../models/nolabs__application__adaptyv_bio__protein_binding_screening__api_models__JobResponse';
import type { nolabs__application__adaptyv_bio__protein_binding_screening__api_models__SetupJobRequest } from '../models/nolabs__application__adaptyv_bio__protein_binding_screening__api_models__SetupJobRequest';
import type { TargetResponse } from '../models/TargetResponse';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class ProteinBindingScreeningService {
    /**
     * Run job
     * @param jobId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static startJobApiV1AdaptyvBioProteinBindingScreeningJobsRunJobIdPost(
        jobId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/adaptyv-bio/protein-binding-screening/jobs/run/{job_id}',
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
     * @returns nolabs__application__adaptyv_bio__protein_binding_screening__api_models__JobResponse Successful Response
     * @throws ApiError
     */
    public static getJobApiV1AdaptyvBioProteinBindingScreeningJobsJobIdGet(
        jobId: string,
    ): CancelablePromise<nolabs__application__adaptyv_bio__protein_binding_screening__api_models__JobResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/adaptyv-bio/protein-binding-screening/jobs/{job_id}',
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
     * @returns nolabs__application__adaptyv_bio__protein_binding_screening__api_models__JobResponse Successful Response
     * @throws ApiError
     */
    public static setupJobApiV1AdaptyvBioProteinBindingScreeningJobsPost(
        requestBody: nolabs__application__adaptyv_bio__protein_binding_screening__api_models__SetupJobRequest,
    ): CancelablePromise<nolabs__application__adaptyv_bio__protein_binding_screening__api_models__JobResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/adaptyv-bio/protein-binding-screening/jobs',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * List targets for the experiment
     * @param searchQuery
     * @returns TargetResponse Successful Response
     * @throws ApiError
     */
    public static listTargetsApiV1AdaptyvBioProteinBindingScreeningListTargetsSearchQueryGet(
        searchQuery: string,
    ): CancelablePromise<Array<TargetResponse>> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/adaptyv-bio/protein-binding-screening/list-targets/{search_query}',
            path: {
                'search_query': searchQuery,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get experiment estimates
     * @param jobId
     * @returns EstimatesResponse Successful Response
     * @throws ApiError
     */
    public static getEstimatesApiV1AdaptyvBioProteinBindingScreeningJobsJobIdEstimatesGet(
        jobId: string,
    ): CancelablePromise<EstimatesResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/adaptyv-bio/protein-binding-screening/jobs/{job_id}/estimates',
            path: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
}
