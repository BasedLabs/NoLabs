/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { ComponentSchema } from './ComponentSchema';
import type { ComponentTemplate_Input } from './ComponentTemplate_Input';
export type WorkflowSchema_Input = {
    workflow_id: string;
    component_templates: Array<ComponentTemplate_Input>;
    components: Array<ComponentSchema>;
    error?: (string | null);
    valid?: boolean;
};

