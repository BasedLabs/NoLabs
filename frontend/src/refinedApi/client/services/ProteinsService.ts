/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { Body_update_protein_api_v1_objects_proteins_patch } from '../models/Body_update_protein_api_v1_objects_proteins_patch';
import type { Body_upload_protein_api_v1_objects_proteins_post } from '../models/Body_upload_protein_api_v1_objects_proteins_post';
import type { ProteinResponse } from '../models/ProteinResponse';
import type { ProteinSearchQuery } from '../models/ProteinSearchQuery';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class ProteinsService {
    /**
     * Search proteins
     * @param requestBody
     * @returns ProteinResponse Successful Response
     * @throws ApiError
     */
    public static searchProteinsApiV1ObjectsProteinsSearchPost(
        requestBody: ProteinSearchQuery,
    ): CancelablePromise<Array<ProteinResponse>> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/objects/proteins/search',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get protein by id
     * @param proteinId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static getProteinApiV1ObjectsProteinsProteinIdGet(
        proteinId: string,
    ): CancelablePromise<(ProteinResponse | null)> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/objects/proteins/{protein_id}',
            path: {
                'protein_id': proteinId,
            },
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
    public static deleteProteinApiV1ObjectsProteinsProteinIdDelete(
        proteinId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'DELETE',
            url: '/api/v1/objects/proteins/{protein_id}',
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
     * @returns ProteinResponse Successful Response
     * @throws ApiError
     */
    public static uploadProteinApiV1ObjectsProteinsPost(
        formData: Body_upload_protein_api_v1_objects_proteins_post,
    ): CancelablePromise<ProteinResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/objects/proteins',
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
     * @returns ProteinResponse Successful Response
     * @throws ApiError
     */
    public static updateProteinApiV1ObjectsProteinsPatch(
        formData: Body_update_protein_api_v1_objects_proteins_patch,
    ): CancelablePromise<ProteinResponse> {
        return __request(OpenAPI, {
            method: 'PATCH',
            url: '/api/v1/objects/proteins',
            formData: formData,
            mediaType: 'multipart/form-data',
            errors: {
                422: `Validation Error`,
            },
        });
    }
}
