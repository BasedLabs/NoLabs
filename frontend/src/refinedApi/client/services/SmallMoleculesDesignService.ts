/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { LogsResponse } from '../models/LogsResponse';
import type { nolabs__application__use_cases__small_molecules_design__api_models__GetJobStatusResponse } from '../models/nolabs__application__use_cases__small_molecules_design__api_models__GetJobStatusResponse';
import type { nolabs__application__use_cases__small_molecules_design__api_models__JobResponse } from '../models/nolabs__application__use_cases__small_molecules_design__api_models__JobResponse';
import type { nolabs__application__use_cases__small_molecules_design__api_models__SetupJobRequest } from '../models/nolabs__application__use_cases__small_molecules_design__api_models__SetupJobRequest';
import type { SmilesResponse } from '../models/SmilesResponse';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class SmallMoleculesDesignService {
    /**
     * Get Job Status
     * @param jobId
     * @returns nolabs__application__use_cases__small_molecules_design__api_models__GetJobStatusResponse Successful Response
     * @throws ApiError
     */
    public static getJobStatusApiV1SmallMoleculesDesignJobsJobIdStatusGet(
        jobId: string,
    ): CancelablePromise<nolabs__application__use_cases__small_molecules_design__api_models__GetJobStatusResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/small-molecules-design/jobs/{job_id}/status',
            path: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Run Learning Stage Job
     * @param jobId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static runLearningStageJobApiV1SmallMoleculesDesignJobsJobIdRunLearningPost(
        jobId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/small-molecules-design/jobs/{job_id}/run/learning',
            path: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Run Sampling Stage Job
     * @param jobId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static runSamplingStageJobApiV1SmallMoleculesDesignJobsJobIdRunSamplingPost(
        jobId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/small-molecules-design/jobs/{job_id}/run/sampling',
            path: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get Job
     * @param jobId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static getJobApiV1SmallMoleculesDesignJobsJobIdGet(
        jobId: string,
    ): CancelablePromise<(nolabs__application__use_cases__small_molecules_design__api_models__JobResponse | null)> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/small-molecules-design/jobs/{job_id}',
            path: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Delete Job
     * @param jobId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static deleteJobApiV1SmallMoleculesDesignJobsJobIdDelete(
        jobId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'DELETE',
            url: '/api/v1/small-molecules-design/jobs/{job_id}',
            path: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Setup Job
     * @param requestBody
     * @returns nolabs__application__use_cases__small_molecules_design__api_models__JobResponse Successful Response
     * @throws ApiError
     */
    public static setupJobApiV1SmallMoleculesDesignJobsPost(
        requestBody: nolabs__application__use_cases__small_molecules_design__api_models__SetupJobRequest,
    ): CancelablePromise<nolabs__application__use_cases__small_molecules_design__api_models__JobResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/small-molecules-design/jobs',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Stop Job
     * @param jobId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static stopJobApiV1SmallMoleculesDesignJobsJobIdStopPost(
        jobId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/small-molecules-design/jobs/{job_id}/stop',
            path: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get Jog Logs
     * @param jobId
     * @returns LogsResponse Successful Response
     * @throws ApiError
     */
    public static getJogLogsApiV1SmallMoleculesDesignJobsJobIdLogsGet(
        jobId: string,
    ): CancelablePromise<LogsResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/small-molecules-design/jobs/{job_id}/logs',
            path: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get Job Smiles
     * @param jobId
     * @returns SmilesResponse Successful Response
     * @throws ApiError
     */
    public static getJobSmilesApiV1SmallMoleculesDesignJobsJobIdSmilesGet(
        jobId: string,
    ): CancelablePromise<Array<SmilesResponse>> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/small-molecules-design/jobs/{job_id}/smiles',
            path: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
}
