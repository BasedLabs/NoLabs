/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { application__folding__api_models__JobResult } from './application__folding__api_models__JobResult';
import type { FoldingBackendEnum } from './FoldingBackendEnum';
export type application__folding__api_models__JobResponse = {
    job_id: string;
    job_name: string;
    backend: FoldingBackendEnum;
    protein_ids: Array<string>;
    result: Array<application__folding__api_models__JobResult>;
    experiment_id: string;
};

