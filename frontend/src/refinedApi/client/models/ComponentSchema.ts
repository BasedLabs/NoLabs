/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { DefaultSchema } from './DefaultSchema';
import type { MappingSchema } from './MappingSchema';
export type ComponentSchema = {
    name: string;
    component_id: string;
    error?: (string | null);
    mappings?: Array<MappingSchema>;
    defaults?: Array<DefaultSchema>;
    'x'?: number;
    'y'?: number;
};

