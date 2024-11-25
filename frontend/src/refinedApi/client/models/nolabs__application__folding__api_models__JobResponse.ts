/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { FoldingBackendEnum } from './FoldingBackendEnum';
import type { nolabs__application__folding__api_models__JobResult } from './nolabs__application__folding__api_models__JobResult';
export type nolabs__application__folding__api_models__JobResponse = {
    job_id: string;
    job_name: string;
    backend: FoldingBackendEnum;
    protein_id: string;
    result: nolabs__application__folding__api_models__JobResult;
    experiment_id: string;
};

