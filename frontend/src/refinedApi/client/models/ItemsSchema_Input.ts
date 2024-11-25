/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { PropertySchema_Input } from './PropertySchema_Input';
export type ItemsSchema_Input = {
    type?: (string | Array<string> | null);
    properties?: (Record<string, PropertySchema_Input> | null);
    required?: Array<string>;
    description?: (string | null);
    enum?: Array<any>;
    $ref?: any;
    const?: null;
    format?: (string | null);
    default?: null;
    example?: null;
    items?: (ItemsSchema_Input | Array<ItemsSchema_Input> | null);
};

