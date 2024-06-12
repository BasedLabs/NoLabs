/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
export type nolabs__refined__application__use_cases__small_molecules_design__api_models__SetupJobRequest = {
    experiment_id: string;
    protein_id: string;
    center_x: number;
    center_y: number;
    center_z: number;
    size_x: number;
    size_y: number;
    size_z: number;
    batch_size?: number;
    minscore?: number;
    epochs?: number;
    job_id?: (string | null);
    job_name?: (string | null);
    sampling_size?: number;
};

