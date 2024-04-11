/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { IntegratorsRequest } from './IntegratorsRequest';
export type nolabs__api_models__conformations__ExperimentPropertiesResponse = {
    pdb_file: string;
    pdb_file_name: string;
    total_frames: number;
    temperature_k: number;
    take_frame_every: number;
    step_size: number;
    replace_non_standard_residues: boolean;
    add_missing_atoms: boolean;
    add_missing_hydrogens: boolean;
    friction_coeff: number;
    ignore_missing_atoms: boolean;
    integrator?: IntegratorsRequest;
};

