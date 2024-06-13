/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { IntegratorsRequest } from './IntegratorsRequest';
export type nolabs__refined__application__use_cases__conformations__api_models__SetupJobRequest = {
    protein_id: string;
    job_id?: (string | null);
    job_name?: (string | null);
    total_frames?: number;
    temperature_k?: number;
    take_frame_every?: number;
    step_size?: number;
    replace_non_standard_residues?: boolean;
    add_missing_atoms?: boolean;
    add_missing_hydrogens?: boolean;
    friction_coeff?: number;
    ignore_missing_atoms?: boolean;
    integrator?: IntegratorsRequest;
    experiment_id?: (string | null);
};

