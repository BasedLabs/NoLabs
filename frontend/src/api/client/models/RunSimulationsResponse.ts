/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { TimelineResponse } from './TimelineResponse';
export type RunSimulationsResponse = {
    experiment_id: string;
    experiment_name: string;
    timeline: Array<TimelineResponse>;
    pdb_content?: (string | null);
};

