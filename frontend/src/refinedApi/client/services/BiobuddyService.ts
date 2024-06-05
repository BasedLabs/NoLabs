/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { CheckBioBuddyEnabledResponse } from '../models/CheckBioBuddyEnabledResponse';
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
}
