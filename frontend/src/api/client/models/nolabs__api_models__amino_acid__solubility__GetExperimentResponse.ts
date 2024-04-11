/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { nolabs__api_models__amino_acid__solubility__AminoAcidResponse } from './nolabs__api_models__amino_acid__solubility__AminoAcidResponse';
import type { nolabs__api_models__amino_acid__solubility__ExperimentPropertiesResponse } from './nolabs__api_models__amino_acid__solubility__ExperimentPropertiesResponse';
export type nolabs__api_models__amino_acid__solubility__GetExperimentResponse = {
    experiment_id: string;
    experiment_name: string;
    amino_acids: Array<nolabs__api_models__amino_acid__solubility__AminoAcidResponse>;
    properties: nolabs__api_models__amino_acid__solubility__ExperimentPropertiesResponse;
};

