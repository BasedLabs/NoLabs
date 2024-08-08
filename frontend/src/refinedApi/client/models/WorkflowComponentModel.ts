/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { DefaultWorkflowComponentModelValue } from './DefaultWorkflowComponentModelValue';
import type { MappingModel } from './MappingModel';
export type WorkflowComponentModel = {
    name: string;
    component_id: string;
    mappings?: Array<MappingModel>;
    error?: (string | null);
    defaults?: Array<DefaultWorkflowComponentModelValue>;
    'x'?: number;
    'y'?: number;
};

