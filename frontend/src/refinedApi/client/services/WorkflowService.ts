/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { GetComponentResponse } from '../models/GetComponentResponse';
import type { GetJobState } from '../models/GetJobState';
import type { WorkflowSchema_Input } from '../models/WorkflowSchema_Input';
import type { WorkflowSchema_Output } from '../models/WorkflowSchema_Output';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class WorkflowService {
    /**
     * Create workflow schema
     * @param experimentId
     * @returns WorkflowSchema_Output Successful Response
     * @throws ApiError
     */
    public static createWorkflowSchemaApiV1WorkflowExperimentIdPost(
        experimentId: string,
    ): CancelablePromise<WorkflowSchema_Output> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/workflow/{experiment_id}',
            path: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get workflow schema
     * @param experimentId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static getSchemaApiV1WorkflowExperimentIdGet(
        experimentId: string,
    ): CancelablePromise<(WorkflowSchema_Output | null)> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/workflow/{experiment_id}',
            path: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Update workflow schema
     * @param requestBody
     * @returns WorkflowSchema_Output Successful Response
     * @throws ApiError
     */
    public static updateWorkflowSchemaApiV1WorkflowPut(
        requestBody: WorkflowSchema_Input,
    ): CancelablePromise<WorkflowSchema_Output> {
        return __request(OpenAPI, {
            method: 'PUT',
            url: '/api/v1/workflow/',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Start workflow
     * @param experimentId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static startWorkflowApiV1WorkflowExperimentIdStartPost(
        experimentId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/workflow/{experiment_id}/start',
            path: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Start workflow component
     * @param experimentId
     * @param componentId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static startComponentApiV1WorkflowExperimentIdStartComponentIdPost(
        experimentId: string,
        componentId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/workflow/{experiment_id}/start/{component_id}',
            path: {
                'experiment_id': experimentId,
                'component_id': componentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get job state
     * @param jobId
     * @returns GetJobState Successful Response
     * @throws ApiError
     */
    public static getJobStateApiV1WorkflowJobJobIdStateGet(
        jobId: string,
    ): CancelablePromise<GetJobState> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/workflow/job/{job_id}/state',
            path: {
                'job_id': jobId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get state
     * @param componentId
     * @returns GetComponentResponse Successful Response
     * @throws ApiError
     */
    public static getComponentStateApiV1WorkflowComponentComponentIdStateGet(
        componentId: string,
    ): CancelablePromise<GetComponentResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/workflow/component/{component_id}/state',
            path: {
                'component_id': componentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
}
