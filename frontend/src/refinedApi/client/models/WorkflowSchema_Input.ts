/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { ComponentSchema } from './ComponentSchema';
import type { ComponentSchemaTemplate_Input } from './ComponentSchemaTemplate_Input';
export type WorkflowSchema_Input = {
    experiment_id: string;
    component_templates: Array<ComponentSchemaTemplate_Input>;
    components: Array<ComponentSchema>;
    error?: (string | null);
    valid?: boolean;
};

