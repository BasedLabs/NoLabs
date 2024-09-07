/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { ComponentSchema } from './ComponentSchema';
import type { ComponentSchemaTemplate_Output } from './ComponentSchemaTemplate_Output';
export type WorkflowSchema_Output = {
    workflow_id: string;
    component_templates: Array<ComponentSchemaTemplate_Output>;
    components: Array<ComponentSchema>;
    error?: (string | null);
    valid?: boolean;
};

