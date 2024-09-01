/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { ItemsTemplate_Input } from './ItemsTemplate_Input';
export type PropertyTemplate_Input = {
    type?: (string | Array<string> | null);
    properties?: (Record<string, PropertyTemplate_Input> | null);
    required?: Array<string>;
    description?: (string | null);
    enum?: Array<any>;
    const?: null;
    format?: (string | null);
    default?: null;
    example?: null;
    title?: (string | null);
    anyOf?: Array<(PropertyTemplate_Input | Record<string, any>)>;
    ref?: (string | null);
    items?: (ItemsTemplate_Input | Array<ItemsTemplate_Input> | null);
};

