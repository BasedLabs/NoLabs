/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { FoldingBackendEnum } from './FoldingBackendEnum';
import type { nolabs__refined__application__use_cases__folding__api_models__JobResult } from './nolabs__refined__application__use_cases__folding__api_models__JobResult';
export type nolabs__refined__application__use_cases__folding__api_models__JobResponse = {
    job_id: string;
    job_name: string;
    backend: FoldingBackendEnum;
    proteins: Array<string>;
    result: Array<nolabs__refined__application__use_cases__folding__api_models__JobResult>;
};
