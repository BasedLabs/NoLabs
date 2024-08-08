/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { Body_update_protein_api_v1_proteins_patch } from '../models/Body_update_protein_api_v1_proteins_patch';
import type { Body_upload_protein_api_v1_proteins_post } from '../models/Body_upload_protein_api_v1_proteins_post';
import type { ProteinContentResponse } from '../models/ProteinContentResponse';
import type { ProteinMetadataResponse } from '../models/ProteinMetadataResponse';
import type { ProteinSearchMetadataQuery } from '../models/ProteinSearchMetadataQuery';
import type { ProteinSearchQuery } from '../models/ProteinSearchQuery';
import type { UploadProteinResponse } from '../models/UploadProteinResponse';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class ProteinsService {
    /**
     * Search proteins content
     * @param requestBody
     * @returns ProteinContentResponse Successful Response
     * @throws ApiError
     */
    public static searchProteinsApiV1ProteinsSearchContentPost(
        requestBody: ProteinSearchQuery,
    ): CancelablePromise<Array<ProteinContentResponse>> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/proteins/search/content',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Search proteins metadata
     * @param requestBody
     * @returns ProteinMetadataResponse Successful Response
     * @throws ApiError
     */
    public static searchProteinsApiV1ProteinsSearchMetadataPost(
        requestBody: ProteinSearchMetadataQuery,
    ): CancelablePromise<Array<ProteinMetadataResponse>> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/proteins/search/metadata',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get protein content by id
     * @param proteinId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static getProteinContentApiV1ProteinsProteinIdContentGet(
        proteinId: string,
    ): CancelablePromise<(ProteinContentResponse | null)> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/proteins/{protein_id}/content',
            path: {
                'protein_id': proteinId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get protein metadata by id
     * @param proteinId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static getProteinMetadataApiV1ProteinsProteinIdMetadataGet(
        proteinId: string,
    ): CancelablePromise<(ProteinMetadataResponse | null)> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/proteins/{protein_id}/metadata',
            path: {
                'protein_id': proteinId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Upload protein
     * @param formData
     * @returns UploadProteinResponse Successful Response
     * @throws ApiError
     */
    public static uploadProteinApiV1ProteinsPost(
        formData: Body_upload_protein_api_v1_proteins_post,
    ): CancelablePromise<UploadProteinResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/proteins',
            formData: formData,
            mediaType: 'multipart/form-data',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Update protein
     * @param formData
     * @returns ProteinContentResponse Successful Response
     * @throws ApiError
     */
    public static updateProteinApiV1ProteinsPatch(
        formData: Body_update_protein_api_v1_proteins_patch,
    ): CancelablePromise<ProteinContentResponse> {
        return __request(OpenAPI, {
            method: 'PATCH',
            url: '/api/v1/proteins',
            formData: formData,
            mediaType: 'multipart/form-data',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Delete protein
     * @param proteinId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static deleteProteinApiV1ProteinsProteinIdDelete(
        proteinId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'DELETE',
            url: '/api/v1/proteins/{protein_id}',
            path: {
                'protein_id': proteinId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
}
