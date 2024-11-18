/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
export type nolabs__application__proteinmpnn__api_models__SetupJobRequest = {
    experiment_id: string;
    protein_id: string;
    num_seq_per_target?: (number | null);
    sampling_temp?: (number | null);
    seed?: (number | null);
    batch_size?: (number | null);
    is_homomer?: (boolean | null);
    chains_to_design?: (Array<string> | null);
    fixed_positions?: (Record<string, Array<number>> | null);
    job_id?: (string | null);
    job_name?: (string | null);
};

