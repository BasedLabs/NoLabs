/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { ComponentModel_Input } from './ComponentModel_Input';
import type { WorkflowComponentModel } from './WorkflowComponentModel';
export type WorkflowSchemaModel_Input = {
    workflow_id: string;
    components: Array<ComponentModel_Input>;
    workflow_components: Array<WorkflowComponentModel>;
    error?: (string | null);
    valid?: boolean;
};

