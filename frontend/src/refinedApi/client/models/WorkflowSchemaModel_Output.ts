/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { ComponentModel_Output } from './ComponentModel_Output';
import type { WorkflowComponentModel } from './WorkflowComponentModel';
export type WorkflowSchemaModel_Output = {
    workflow_id: string;
    components: Array<ComponentModel_Output>;
    workflow_components: Array<WorkflowComponentModel>;
    error?: (string | null);
    valid?: boolean;
};

