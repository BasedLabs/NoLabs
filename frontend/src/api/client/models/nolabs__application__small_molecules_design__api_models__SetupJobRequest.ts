/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
export type nolabs__application__small_molecules_design__api_models__SetupJobRequest = {
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
    experiment_id?: (string | null);
};

