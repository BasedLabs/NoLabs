/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { nolabs__api_models__conformations__ExperimentPropertiesResponse } from './nolabs__api_models__conformations__ExperimentPropertiesResponse';
import type { TimelineResponse } from './TimelineResponse';
export type nolabs__api_models__conformations__GetExperimentResponse = {
    experiment_id: string;
    experiment_name: string;
    properties: nolabs__api_models__conformations__ExperimentPropertiesResponse;
    timeline: Array<TimelineResponse>;
    pdb_file: (string | null);
};

