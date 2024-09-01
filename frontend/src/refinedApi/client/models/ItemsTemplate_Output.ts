/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { PropertyTemplate_Output } from './PropertyTemplate_Output';
export type ItemsTemplate_Output = {
    type?: (string | Array<string> | null);
    properties?: (Record<string, PropertyTemplate_Output> | null);
    required?: Array<string>;
    description?: (string | null);
    enum?: Array<any>;
    $ref?: any;
    const?: null;
    format?: (string | null);
    default?: null;
    example?: null;
    items?: (ItemsTemplate_Output | Array<ItemsTemplate_Output> | null);
};

