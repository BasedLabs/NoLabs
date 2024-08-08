/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { FoldingBackendEnum } from './FoldingBackendEnum';
export type nolabs__application__use_cases__folding__api_models__SetupJobRequest = {
    experiment_id: string;
    backend: (FoldingBackendEnum | null);
    protein_ids: Array<string>;
    job_id?: (string | null);
    job_name?: (string | null);
};

