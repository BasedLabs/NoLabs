/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { ItemsSchema_Input } from './ItemsSchema_Input';
export type PropertySchema_Input = {
    type?: (string | Array<string> | null);
    properties?: (Record<string, PropertySchema_Input> | null);
    required?: Array<string>;
    description?: (string | null);
    enum?: Array<any>;
    const?: null;
    format?: (string | null);
    default?: null;
    example?: null;
    title?: (string | null);
    anyOf?: Array<(PropertySchema_Input | Record<string, any>)>;
    ref?: (string | null);
    items?: (ItemsSchema_Input | Array<ItemsSchema_Input> | null);
};

