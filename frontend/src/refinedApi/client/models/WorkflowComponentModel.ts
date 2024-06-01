/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { DefaultWorkflowComponentModelValue } from './DefaultWorkflowComponentModelValue';
import type { JobValidationError } from './JobValidationError';
import type { MappingModel } from './MappingModel';
export type WorkflowComponentModel = {
    name: string;
    component_id: string;
    job_ids?: Array<string>;
    mappings?: Array<MappingModel>;
    error?: (string | null);
    defaults?: Array<DefaultWorkflowComponentModelValue>;
    jobs_errors: Array<JobValidationError>;
};

