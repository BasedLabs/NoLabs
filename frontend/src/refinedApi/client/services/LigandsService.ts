/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { Body_update_ligand_api_v1_objects_ligands_patch } from '../models/Body_update_ligand_api_v1_objects_ligands_patch';
import type { Body_upload_ligand_api_v1_objects_ligands_post } from '../models/Body_upload_ligand_api_v1_objects_ligands_post';
import type { LigandContentResponse } from '../models/LigandContentResponse';
import type { LigandMetadataResponse } from '../models/LigandMetadataResponse';
import type { LigandSearchContentQuery } from '../models/LigandSearchContentQuery';
import type { LigandSearchMetadataQuery } from '../models/LigandSearchMetadataQuery';
import type { UploadLigandResponse } from '../models/UploadLigandResponse';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class LigandsService {
    /**
     * Search ligands content
     * @param requestBody
     * @returns LigandContentResponse Successful Response
     * @throws ApiError
     */
    public static searchLigandsApiV1ObjectsLigandsSearchContentPost(
        requestBody: LigandSearchContentQuery,
    ): CancelablePromise<Array<LigandContentResponse>> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/objects/ligands/search/content',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Search ligands metadata
     * @param requestBody
     * @returns LigandMetadataResponse Successful Response
     * @throws ApiError
     */
    public static searchLigandsApiV1ObjectsLigandsSearchMetadataPost(
        requestBody: LigandSearchMetadataQuery,
    ): CancelablePromise<Array<LigandMetadataResponse>> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/objects/ligands/search/metadata',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get ligand content by id
     * @param ligandId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static getLigandContentApiV1ObjectsLigandsLigandIdContentGet(
        ligandId: string,
    ): CancelablePromise<(LigandContentResponse | null)> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/objects/ligands/{ligand_id}/content',
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
     * @returns UploadLigandResponse Successful Response
     * @throws ApiError
     */
    public static uploadLigandApiV1ObjectsLigandsPost(
        formData: Body_upload_ligand_api_v1_objects_ligands_post,
    ): CancelablePromise<UploadLigandResponse> {
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
     * @returns LigandContentResponse Successful Response
     * @throws ApiError
     */
    public static updateLigandApiV1ObjectsLigandsPatch(
        formData: Body_update_ligand_api_v1_objects_ligands_patch,
    ): CancelablePromise<LigandContentResponse> {
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
}
