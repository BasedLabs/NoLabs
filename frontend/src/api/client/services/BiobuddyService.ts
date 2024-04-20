/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { CheckBioBuddyEnabledResponse } from '../models/CheckBioBuddyEnabledResponse';
import type { LoadConversationResponse } from '../models/LoadConversationResponse';
import type { SaveMessageResponse } from '../models/SaveMessageResponse';
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
     * Save Message
     * @param experimentId
     * @param messageContent
     * @returns SaveMessageResponse Successful Response
     * @throws ApiError
     */
    public static saveMessageApiV1BiobuddySaveMessagePost(
        experimentId: string,
        messageContent: string,
    ): CancelablePromise<SaveMessageResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/biobuddy/save-message',
            query: {
                'experiment_id': experimentId,
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
     * @param messageContent
     * @returns SendQueryResponse Successful Response
     * @throws ApiError
     */
    public static sendQueryApiV1BiobuddySendQueryPost(
        experimentId: string,
        messageContent: string,
    ): CancelablePromise<SendQueryResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/biobuddy/send-query',
            query: {
                'experiment_id': experimentId,
                'message_content': messageContent,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
}
