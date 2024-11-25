/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { nolabs__application__proteinmpnn__api_models__JobResult } from './nolabs__application__proteinmpnn__api_models__JobResult';
export type nolabs__application__proteinmpnn__api_models__JobResponse = {
    job_id: string;
    job_name: string;
    num_seq_per_target: number;
    sampling_temp: number;
    seed: number;
    batch_size: number;
    is_homomer: boolean;
    protein_id: string;
    result: Array<nolabs__application__proteinmpnn__api_models__JobResult>;
    chains_to_design?: (Array<string> | null);
    fixed_positions?: (Record<string, Array<number>> | null);
};

