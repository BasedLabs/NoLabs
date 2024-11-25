/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { ItemsSchema_Output } from './ItemsSchema_Output';
export type PropertySchema_Output = {
    type?: (string | Array<string> | null);
    properties?: (Record<string, PropertySchema_Output> | null);
    required?: Array<string>;
    description?: (string | null);
    enum?: Array<any>;
    const?: null;
    format?: (string | null);
    default?: null;
    example?: null;
    title?: (string | null);
    anyOf?: Array<(PropertySchema_Output | Record<string, any>)>;
    ref?: (string | null);
    items?: (ItemsSchema_Output | Array<ItemsSchema_Output> | null);
};

