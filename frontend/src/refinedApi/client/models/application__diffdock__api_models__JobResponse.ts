/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { application__diffdock__api_models__JobResult } from './application__diffdock__api_models__JobResult';
export type application__diffdock__api_models__JobResponse = {
    job_id: string;
    job_name: string;
    samples_per_complex: number;
    protein_id: string;
    ligand_id: string;
    result: Array<application__diffdock__api_models__JobResult>;
};

