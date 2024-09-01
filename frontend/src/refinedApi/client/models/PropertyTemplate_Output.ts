/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { ItemsTemplate_Output } from './ItemsTemplate_Output';
export type PropertyTemplate_Output = {
    type?: (string | Array<string> | null);
    properties?: (Record<string, PropertyTemplate_Output> | null);
    required?: Array<string>;
    description?: (string | null);
    enum?: Array<any>;
    const?: null;
    format?: (string | null);
    default?: null;
    example?: null;
    title?: (string | null);
    anyOf?: Array<(PropertyTemplate_Output | Record<string, any>)>;
    ref?: (string | null);
    items?: (ItemsTemplate_Output | Array<ItemsTemplate_Output> | null);
};

