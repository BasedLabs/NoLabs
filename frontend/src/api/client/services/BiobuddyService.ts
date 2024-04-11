/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { CheckBioBuddyEnabledResponse } from '../models/CheckBioBuddyEnabledResponse';
import type { LoadConversationResponse } from '../models/LoadConversationResponse';
import type { SendMessageResponse } from '../models/SendMessageResponse';
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
     * Send Message
     * @param experimentId
     * @param messageContent
     * @returns SendMessageResponse Successful Response
     * @throws ApiError
     */
    public static sendMessageApiV1BiobuddySendMessagePost(
        experimentId: string,
        messageContent: string,
    ): CancelablePromise<SendMessageResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/biobuddy/send-message',
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
