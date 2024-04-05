/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { nolabs__api_models__small_molecules_design__ExperimentPropertiesResponse } from './nolabs__api_models__small_molecules_design__ExperimentPropertiesResponse';
export type nolabs__api_models__small_molecules_design__GetExperimentResponse = {
    experiment_id: string;
    experiment_name: string;
    created_at: string;
    running: boolean;
    learning_completed: boolean;
    properties: nolabs__api_models__small_molecules_design__ExperimentPropertiesResponse;
};

