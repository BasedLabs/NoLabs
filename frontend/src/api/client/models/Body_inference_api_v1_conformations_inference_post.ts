/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { IntegratorsRequest } from './IntegratorsRequest';
export type Body_inference_api_v1_conformations_inference_post = {
    pdb_file: Blob;
    experiment_name: string;
    experiment_id?: string;
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
};

