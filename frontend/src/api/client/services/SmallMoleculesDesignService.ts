/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { Body_save_properties_api_v1_small_molecules_design_experiment__experiment_id__props_post } from '../models/Body_save_properties_api_v1_small_molecules_design_experiment__experiment_id__props_post';
import type { ChangeExperimentNameRequest } from '../models/ChangeExperimentNameRequest';
import type { ExperimentMetadataResponse } from '../models/ExperimentMetadataResponse';
import type { GetExperimentStatusResponse } from '../models/GetExperimentStatusResponse';
import type { LogsResponse } from '../models/LogsResponse';
import type { nolabs__api_models__small_molecules_design__GetExperimentResponse } from '../models/nolabs__api_models__small_molecules_design__GetExperimentResponse';
import type { SamplingSizeRequest } from '../models/SamplingSizeRequest';
import type { SmilesResponse } from '../models/SmilesResponse';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class SmallMoleculesDesignService {
    /**
     * Status
     * @param experimentId
     * @returns GetExperimentStatusResponse Successful Response
     * @throws ApiError
     */
    public static statusApiV1SmallMoleculesDesignExperimentExperimentIdStatusGet(
        experimentId: string,
    ): CancelablePromise<GetExperimentStatusResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/small-molecules-design/experiment/{experiment_id}/status',
            path: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Learning
     * @param experimentId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static learningApiV1SmallMoleculesDesignExperimentExperimentIdLearningPost(
        experimentId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/small-molecules-design/experiment/{experiment_id}/learning',
            path: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Sampling
     * @param experimentId
     * @param requestBody
     * @returns any Successful Response
     * @throws ApiError
     */
    public static samplingApiV1SmallMoleculesDesignExperimentExperimentIdSamplingPost(
        experimentId: string,
        requestBody: SamplingSizeRequest,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/small-molecules-design/experiment/{experiment_id}/sampling',
            path: {
                'experiment_id': experimentId,
            },
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Get Experiment
     * @param experimentId
     * @returns nolabs__api_models__small_molecules_design__GetExperimentResponse Successful Response
     * @throws ApiError
     */
    public static getExperimentApiV1SmallMoleculesDesignExperimentExperimentIdGet(
        experimentId: string,
    ): CancelablePromise<nolabs__api_models__small_molecules_design__GetExperimentResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/small-molecules-design/experiment/{experiment_id}',
            path: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Delete
     * @param experimentId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static deleteApiV1SmallMoleculesDesignExperimentExperimentIdDelete(
        experimentId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'DELETE',
            url: '/api/v1/small-molecules-design/experiment/{experiment_id}',
            path: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Save Properties
     * @param experimentId
     * @param centerX
     * @param centerY
     * @param centerZ
     * @param sizeX
     * @param sizeY
     * @param sizeZ
     * @param formData
     * @param epochs
     * @param batchSize
     * @param minscore
     * @returns any Successful Response
     * @throws ApiError
     */
    public static savePropertiesApiV1SmallMoleculesDesignExperimentExperimentIdPropsPost(
        experimentId: string,
        centerX: number,
        centerY: number,
        centerZ: number,
        sizeX: number,
        sizeY: number,
        sizeZ: number,
        formData: Body_save_properties_api_v1_small_molecules_design_experiment__experiment_id__props_post,
        epochs: number = 50,
        batchSize: number = 128,
        minscore: number = 0.4,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/small-molecules-design/experiment/{experiment_id}/props',
            path: {
                'experiment_id': experimentId,
            },
            query: {
                'center_x': centerX,
                'center_y': centerY,
                'center_z': centerZ,
                'size_x': sizeX,
                'size_y': sizeY,
                'size_z': sizeZ,
                'epochs': epochs,
                'batch_size': batchSize,
                'minscore': minscore,
            },
            formData: formData,
            mediaType: 'multipart/form-data',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Stop
     * @param experimentId
     * @returns any Successful Response
     * @throws ApiError
     */
    public static stopApiV1SmallMoleculesDesignExperimentExperimentIdStopPost(
        experimentId: string,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/small-molecules-design/experiment/{experiment_id}/stop',
            path: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Logs
     * @param experimentId
     * @returns LogsResponse Successful Response
     * @throws ApiError
     */
    public static logsApiV1SmallMoleculesDesignExperimentExperimentIdLogsGet(
        experimentId: string,
    ): CancelablePromise<LogsResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/small-molecules-design/experiment/{experiment_id}/logs',
            path: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Smiles
     * @param experimentId
     * @returns SmilesResponse Successful Response
     * @throws ApiError
     */
    public static smilesApiV1SmallMoleculesDesignExperimentExperimentIdSmilesGet(
        experimentId: string,
    ): CancelablePromise<Array<SmilesResponse>> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/small-molecules-design/experiment/{experiment_id}/smiles',
            path: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Experiments
     * @returns ExperimentMetadataResponse Successful Response
     * @throws ApiError
     */
    public static experimentsApiV1SmallMoleculesDesignExperimentsGet(): CancelablePromise<Array<ExperimentMetadataResponse>> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/small-molecules-design/experiments',
        });
    }
    /**
     * Change Experiment Name
     * @param requestBody
     * @returns any Successful Response
     * @throws ApiError
     */
    public static changeExperimentNameApiV1SmallMoleculesDesignExperimentNamePost(
        requestBody: ChangeExperimentNameRequest,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/small-molecules-design/experiment/name',
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Create Experiment
     * @returns ExperimentMetadataResponse Successful Response
     * @throws ApiError
     */
    public static createExperimentApiV1SmallMoleculesDesignExperimentCreatePost(): CancelablePromise<ExperimentMetadataResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/small-molecules-design/experiment/create',
        });
    }
}
