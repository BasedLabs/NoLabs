/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { CheckBioBuddyEnabledResponse } from '../models/CheckBioBuddyEnabledResponse';
import type { CreateFunctionCallMessageResponse } from '../models/CreateFunctionCallMessageResponse';
import type { CreateMessageResponse } from '../models/CreateMessageResponse';
import type { EditMessageResponse } from '../models/EditMessageResponse';
import type { FunctionCall_Input } from '../models/FunctionCall_Input';
import type { GetAvailableFunctionCallsResponse } from '../models/GetAvailableFunctionCallsResponse';
import type { LoadConversationResponse } from '../models/LoadConversationResponse';
import type { SendQueryResponse } from '../models/SendQueryResponse';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class BiobuddyService {
    /**
     * Check Biobuddy Enabled
     * @returns CheckBioBuddyEnabledResponse Successful Response
     * @throws ApiError
     */
    public static checkBiobuddyEnabledApiV1BiobuddyCheckBiobuddyEnabledGet(): CancelablePromise<CheckBioBuddyEnabledResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/biobuddy/check_biobuddy_enabled',
        });
    }
    /**
     * Load Conversation
     * @param experimentId
     * @returns LoadConversationResponse Successful Response
     * @throws ApiError
     */
    public static loadConversationApiV1BiobuddyLoadConversationGet(
        experimentId: string,
    ): CancelablePromise<LoadConversationResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/biobuddy/load-conversation',
            query: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Create Message
     * @param experimentId
     * @param messageId
     * @param messageContent
     * @param role
     * @returns CreateMessageResponse Successful Response
     * @throws ApiError
     */
    public static createMessageApiV1BiobuddyMessageCreatePost(
        experimentId: string,
        messageId: string,
        messageContent: string,
        role: string,
    ): CancelablePromise<CreateMessageResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/biobuddy/message/create',
            query: {
                'experiment_id': experimentId,
                'message_id': messageId,
                'message_content': messageContent,
                'role': role,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Create Function Call Message
     * @param experimentId
     * @param messageId
     * @param role
     * @param requestBody
     * @returns CreateFunctionCallMessageResponse Successful Response
     * @throws ApiError
     */
    public static createFunctionCallMessageApiV1BiobuddyFunctionCreatePost(
        experimentId: string,
        messageId: string,
        role: string,
        requestBody: FunctionCall_Input,
    ): CancelablePromise<CreateFunctionCallMessageResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/biobuddy/function/create',
            query: {
                'experiment_id': experimentId,
                'message_id': messageId,
                'role': role,
            },
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Edit Message
     * @param experimentId
     * @param messageId
     * @param messageContent
     * @returns EditMessageResponse Successful Response
     * @throws ApiError
     */
    public static editMessageApiV1BiobuddyMessageEditPost(
        experimentId: string,
        messageId: string,
        messageContent: string,
    ): CancelablePromise<EditMessageResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/biobuddy/message/edit',
            query: {
                'experiment_id': experimentId,
                'message_id': messageId,
                'message_content': messageContent,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Send Query
     * @param experimentId
     * @param actionText
     * @returns SendQueryResponse Successful Response
     * @throws ApiError
     */
    public static sendQueryApiV1BiobuddyQueryPost(
        experimentId: string,
        actionText: string,
    ): CancelablePromise<SendQueryResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/biobuddy/query',
            query: {
                'experiment_id': experimentId,
                'action_text': actionText,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get Tools
     * @returns GetAvailableFunctionCallsResponse Successful Response
     * @throws ApiError
     */
    public static getToolsApiV1BiobuddyToolsGet(): CancelablePromise<GetAvailableFunctionCallsResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/biobuddy/tools',
        });
    }
}
