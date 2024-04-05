/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { Body_inference_api_v1_gene_ontology_inference_post } from '../models/Body_inference_api_v1_gene_ontology_inference_post';
import type { ChangeExperimentNameRequest } from '../models/ChangeExperimentNameRequest';
import type { ExperimentMetadataResponse } from '../models/ExperimentMetadataResponse';
import type { nolabs__api_models__amino_acid__gene_ontology__GetExperimentResponse } from '../models/nolabs__api_models__amino_acid__gene_ontology__GetExperimentResponse';
import type { RunGeneOntologyResponse } from '../models/RunGeneOntologyResponse';
import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';
export class GeneOntologyService {
    /**
     * Inference
     * @param formData
     * @returns RunGeneOntologyResponse Successful Response
     * @throws ApiError
     */
    public static inferenceApiV1GeneOntologyInferencePost(
        formData: Body_inference_api_v1_gene_ontology_inference_post,
    ): CancelablePromise<RunGeneOntologyResponse> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/gene-ontology/inference',
            formData: formData,
            mediaType: 'multipart/form-data',
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
    public static experimentsApiV1GeneOntologyExperimentsGet(): CancelablePromise<Array<ExperimentMetadataResponse>> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/gene-ontology/experiments',
        });
    }
    /**
     * Get Experiment
     * @param experimentId
     * @returns nolabs__api_models__amino_acid__gene_ontology__GetExperimentResponse Successful Response
     * @throws ApiError
     */
    public static getExperimentApiV1GeneOntologyGetExperimentGet(
        experimentId: string,
    ): CancelablePromise<nolabs__api_models__amino_acid__gene_ontology__GetExperimentResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/gene-ontology/get-experiment',
            query: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Delete Experiment
     * @param experimentId
     * @returns nolabs__api_models__amino_acid__gene_ontology__GetExperimentResponse Successful Response
     * @throws ApiError
     */
    public static deleteExperimentApiV1GeneOntologyDeleteExperimentDelete(
        experimentId: string,
    ): CancelablePromise<nolabs__api_models__amino_acid__gene_ontology__GetExperimentResponse> {
        return __request(OpenAPI, {
            method: 'DELETE',
            url: '/api/v1/gene-ontology/delete-experiment',
            query: {
                'experiment_id': experimentId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }
    /**
     * Change Experiment Name
     * @param requestBody
     * @returns any Successful Response
     * @throws ApiError
     */
    public static changeExperimentNameApiV1GeneOntologyChangeExperimentNamePost(
        requestBody: ChangeExperimentNameRequest,
    ): CancelablePromise<any> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/api/v1/gene-ontology/change-experiment-name',
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
    public static createExperimentApiV1GeneOntologyCreateExperimentGet(): CancelablePromise<ExperimentMetadataResponse> {
        return __request(OpenAPI, {
            method: 'GET',
            url: '/api/v1/gene-ontology/create-experiment',
        });
    }
}
