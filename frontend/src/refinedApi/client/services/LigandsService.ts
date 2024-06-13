/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { Body_update_ligand_api_v1_objects_ligands_patch } from '../models/Body_update_ligand_api_v1_objects_ligands_patch';
import type { Body_upload_ligand_api_v1_objects_ligands_post } from '../models/Body_upload_ligand_api_v1_objects_ligands_post';
import type { LigandResponse } from '../models/LigandResponse';
import type { LigandSearchQuery } from '../models/LigandSearchQuery';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class LigandsService {
    /**
     * Search ligands
     * @param requestBody
     * @returns LigandResponse Successful Response
     * @throws ApiError
     */
    public static searchLigandsApiV1ObjectsLigandsSearchPost(
        requestBody: LigandSearchQuery,
    ): CancelablePromise<Array<LigandResponse>> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/objects/ligands/search',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get ligand by id
     * @param ligandId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static getLigandApiV1ObjectsLigandsLigandIdGet(
        ligandId: string,
    ): CancelablePromise<(LigandResponse | null)> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/objects/ligands/{ligand_id}',
            path: {
                'ligand_id': ligandId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Delete ligand
     * @param ligandId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static deleteLigandApiV1ObjectsLigandsLigandIdDelete(
        ligandId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'DELETE',
            url: '/api/v1/objects/ligands/{ligand_id}',
            path: {
                'ligand_id': ligandId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Upload ligand
     * @param formData
     * @returns LigandResponse Successful Response
     * @throws ApiError
     */
    public static uploadLigandApiV1ObjectsLigandsPost(
        formData: Body_upload_ligand_api_v1_objects_ligands_post,
    ): CancelablePromise<LigandResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/objects/ligands',
            formData: formData,
            mediaType: 'multipart/form-data',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Update ligand
     * @param formData
     * @returns LigandResponse Successful Response
     * @throws ApiError
     */
    public static updateLigandApiV1ObjectsLigandsPatch(
        formData: Body_update_ligand_api_v1_objects_ligands_patch,
    ): CancelablePromise<LigandResponse> {
        return __request(OpenAPI, {
            method: 'PATCH',
            url: '/api/v1/objects/ligands',
            formData: formData,
            mediaType: 'multipart/form-data',
            errors: {
                422: `Validation Error`,
            },
        });
    }
}
