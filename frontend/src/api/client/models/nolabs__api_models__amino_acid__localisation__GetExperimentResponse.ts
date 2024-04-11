/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { nolabs__api_models__amino_acid__localisation__AminoAcidResponse } from './nolabs__api_models__amino_acid__localisation__AminoAcidResponse';
import type { nolabs__api_models__amino_acid__localisation__ExperimentPropertiesResponse } from './nolabs__api_models__amino_acid__localisation__ExperimentPropertiesResponse';
export type nolabs__api_models__amino_acid__localisation__GetExperimentResponse = {
    experiment_id: string;
    experiment_name: string;
    amino_acids: Array<nolabs__api_models__amino_acid__localisation__AminoAcidResponse>;
    properties: nolabs__api_models__amino_acid__localisation__ExperimentPropertiesResponse;
};

