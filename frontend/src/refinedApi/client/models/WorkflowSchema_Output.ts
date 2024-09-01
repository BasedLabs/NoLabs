/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { ComponentSchema } from './ComponentSchema';
import type { ComponentTemplate_Output } from './ComponentTemplate_Output';
export type WorkflowSchema_Output = {
    workflow_id: string;
    component_templates: Array<ComponentTemplate_Output>;
    components: Array<ComponentSchema>;
    error?: (string | null);
    valid?: boolean;
};

