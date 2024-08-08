/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { nolabs__application__use_cases__blast__api_models__JobResult } from './nolabs__application__use_cases__blast__api_models__JobResult';
export type nolabs__application__use_cases__blast__api_models__JobResponse = {
    job_id: string;
    job_name: string;
    protein_id: string;
    descriptions: number;
    alignments: number;
    hitlist_size: number;
    expect: number;
    result: Array<nolabs__application__use_cases__blast__api_models__JobResult>;
};

