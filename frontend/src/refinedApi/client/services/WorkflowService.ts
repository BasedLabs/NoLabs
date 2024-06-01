/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { AllWorkflowSchemasResponse } from '../models/AllWorkflowSchemasResponse';
import type { WorkflowSchemaModel_Input } from '../models/WorkflowSchemaModel_Input';
import type { WorkflowSchemaModel_Output } from '../models/WorkflowSchemaModel_Output';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class WorkflowService {
    /**
     * Create workflow schema
     * @param experimentId
     * @returns WorkflowSchemaModel_Output Successful Response
     * @throws ApiError
     */
    public static createSchemaApiV1WorkflowExperimentIdPost(
        experimentId: string,
    ): CancelablePromise<WorkflowSchemaModel_Output> {
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
     * Delete workflow schema
     * @param workflowId
     * @returns WorkflowSchemaModel_Output Successful Response
     * @throws ApiError
     */
    public static createSchemaApiV1WorkflowWorkflowIdDelete(
        workflowId: string,
    ): CancelablePromise<WorkflowSchemaModel_Output> {
        return __request(OpenAPI, {
            method: 'DELETE',
            url: '/api/v1/workflow/{workflow_id}',
            path: {
                'workflow_id': workflowId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get workflow schema
     * @param workflowId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static getSchemaApiV1WorkflowWorkflowIdGet(
        workflowId: string,
    ): CancelablePromise<(WorkflowSchemaModel_Output | null)> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/workflow/{workflow_id}',
            path: {
                'workflow_id': workflowId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * All workflow schemas
     * @param experimentId
     * @returns AllWorkflowSchemasResponse Successful Response
     * @throws ApiError
     */
    public static getSchemaApiV1WorkflowAllExperimentIdGet(
        experimentId: string,
    ): CancelablePromise<AllWorkflowSchemasResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/workflow/all/{experiment_id}',
            path: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Set workflow schema
     * @param requestBody
     * @returns WorkflowSchemaModel_Output Successful Response
     * @throws ApiError
     */
    public static setWorkflowSchemaApiV1WorkflowPut(
        requestBody: WorkflowSchemaModel_Input,
    ): CancelablePromise<WorkflowSchemaModel_Output> {
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
     * Start workflow schema
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
}
