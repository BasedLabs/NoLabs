/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { ComponentStateEnum } from './ComponentStateEnum';
import type { PropertyErrorResponse } from './PropertyErrorResponse';
export type GetComponentResponse = {
    id: string;
    input_schema: Record<string, any>;
    output_schema: Record<string, any>;
    input_value_dict: Record<string, any>;
    output_value_dict: Record<string, any>;
    previous_component_ids: Array<string>;
    input_errors: Array<PropertyErrorResponse>;
    output_errors: Array<PropertyErrorResponse>;
    state: ComponentStateEnum;
    job_ids: Array<string>;
    state_message?: (string | null);
};

