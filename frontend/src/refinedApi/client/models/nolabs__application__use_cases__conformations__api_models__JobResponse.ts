/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { IntegratorsRequest } from './IntegratorsRequest';
import type { TimelineResponse } from './TimelineResponse';
export type nolabs__application__use_cases__conformations__api_models__JobResponse = {
    job_id: string;
    job_name: string;
    experiment_id: string;
    protein_id: string;
    timeline: Array<TimelineResponse>;
    total_frames: number;
    temperature_k: number;
    take_frame_every: number;
    step_size: number;
    replace_non_standard_residues: boolean;
    add_missing_atoms: boolean;
    add_missing_hydrogens: boolean;
    friction_coeff: number;
    ignore_missing_atoms: boolean;
    integrator: IntegratorsRequest;
    md_content?: (Blob | null);
};

