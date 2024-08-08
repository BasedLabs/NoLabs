/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { ItemsModel_Input } from './ItemsModel_Input';
export type PropertyModel_Input = {
    type?: (string | Array<string> | null);
    properties?: (Record<string, PropertyModel_Input> | null);
    required?: Array<string>;
    description?: (string | null);
    enum?: Array<any>;
    const?: null;
    format?: (string | null);
    default?: null;
    example?: null;
    title?: (string | null);
    anyOf?: Array<(PropertyModel_Input | Record<string, any>)>;
    ref?: (string | null);
    items?: (ItemsModel_Input | Array<ItemsModel_Input> | null);
};

