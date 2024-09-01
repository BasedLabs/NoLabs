/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { PropertyTemplate_Input } from './PropertyTemplate_Input';
export type ItemsTemplate_Input = {
    type?: (string | Array<string> | null);
    properties?: (Record<string, PropertyTemplate_Input> | null);
    required?: Array<string>;
    description?: (string | null);
    enum?: Array<any>;
    $ref?: any;
    const?: null;
    format?: (string | null);
    default?: null;
    example?: null;
    items?: (ItemsTemplate_Input | Array<ItemsTemplate_Input> | null);
};

