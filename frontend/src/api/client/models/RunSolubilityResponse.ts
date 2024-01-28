/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { nolabs__api_models__solubility__AminoAcidResponse } from './nolabs__api_models__solubility__AminoAcidResponse';
export type RunSolubilityResponse = {
    experiment_id: string;
    experiment_name: string;
    amino_acids: Array<nolabs__api_models__solubility__AminoAcidResponse>;
    errors?: Array<string>;
};

