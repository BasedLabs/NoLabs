/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { InputPropertyErrorResponse } from './InputPropertyErrorResponse';
import type { JobErrorResponse } from './JobErrorResponse';
export type GetComponentStateResponse = {
    input_dict: Record<string, any>;
    output_dict: Record<string, any>;
    job_ids: Array<string>;
    input_property_errors: Array<InputPropertyErrorResponse>;
    last_exceptions: Array<string>;
    jobs_errors?: Array<JobErrorResponse>;
};

