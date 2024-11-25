/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { EstimatesResponse } from '../models/EstimatesResponse';
import type { nolabs__application__adaptyv_bio__protein_affinity_characterization__api_models__JobResponse } from '../models/nolabs__application__adaptyv_bio__protein_affinity_characterization__api_models__JobResponse';
import type { nolabs__application__adaptyv_bio__protein_affinity_characterization__api_models__SetupJobRequest } from '../models/nolabs__application__adaptyv_bio__protein_affinity_characterization__api_models__SetupJobRequest';
import type { TargetResponse } from '../models/TargetResponse';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class ProteinAffinityCharacterizationService {
    /**
     * Run job
     * @param jobId
     * @returns nolabs__application__adaptyv_bio__protein_affinity_characterization__api_models__JobResponse Successful Response
     * @throws ApiError
     */
    public static startJobApiV1AdaptyvBioProteinAffinityJobsRunJobIdPost(
        jobId: string,
    ): CancelablePromise<nolabs__application__adaptyv_bio__protein_affinity_characterization__api_models__JobResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/adaptyv-bio/protein-affinity/jobs/run/{job_id}',
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
     * @returns nolabs__application__adaptyv_bio__protein_affinity_characterization__api_models__JobResponse Successful Response
     * @throws ApiError
     */
    public static getJobApiV1AdaptyvBioProteinAffinityJobsJobIdGet(
        jobId: string,
    ): CancelablePromise<nolabs__application__adaptyv_bio__protein_affinity_characterization__api_models__JobResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/adaptyv-bio/protein-affinity/jobs/{job_id}',
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
     * @returns nolabs__application__adaptyv_bio__protein_affinity_characterization__api_models__JobResponse Successful Response
     * @throws ApiError
     */
    public static setupJobApiV1AdaptyvBioProteinAffinityJobsPost(
        requestBody: nolabs__application__adaptyv_bio__protein_affinity_characterization__api_models__SetupJobRequest,
    ): CancelablePromise<nolabs__application__adaptyv_bio__protein_affinity_characterization__api_models__JobResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/adaptyv-bio/protein-affinity/jobs',
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
    public static listTargetsApiV1AdaptyvBioProteinAffinityListTargetsSearchQueryGet(
        searchQuery: string,
    ): CancelablePromise<Array<TargetResponse>> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/adaptyv-bio/protein-affinity/list-targets/{search_query}',
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
    public static getEstimatesApiV1AdaptyvBioProteinAffinityJobsJobIdEstimatesGet(
        jobId: string,
    ): CancelablePromise<EstimatesResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/adaptyv-bio/protein-affinity/jobs/{job_id}/estimates',
            path: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
}
