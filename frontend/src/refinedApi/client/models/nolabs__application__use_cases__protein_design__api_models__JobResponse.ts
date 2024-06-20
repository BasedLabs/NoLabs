/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
export type nolabs__application__use_cases__protein_design__api_models__JobResponse = {
    job_id: string;
    job_name: string;
    experiment_id: string;
    protein_id: string;
    binder_ids: Array<string>;
    contig: string;
    number_of_designs?: number;
    timesteps?: (number | null);
    hotspots?: (string | null);
};

