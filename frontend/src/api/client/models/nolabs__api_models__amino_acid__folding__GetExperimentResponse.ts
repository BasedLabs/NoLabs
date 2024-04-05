/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { nolabs__api_models__amino_acid__folding__AminoAcidResponse } from './nolabs__api_models__amino_acid__folding__AminoAcidResponse';
import type { nolabs__api_models__amino_acid__folding__ExperimentPropertiesResponse } from './nolabs__api_models__amino_acid__folding__ExperimentPropertiesResponse';
export type nolabs__api_models__amino_acid__folding__GetExperimentResponse = {
    experiment_id: string;
    experiment_name: string;
    amino_acids: Array<nolabs__api_models__amino_acid__folding__AminoAcidResponse>;
    properties: nolabs__api_models__amino_acid__folding__ExperimentPropertiesResponse;
};

