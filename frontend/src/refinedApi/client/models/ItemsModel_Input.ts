/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { PropertyModel_Input } from './PropertyModel_Input';
export type ItemsModel_Input = {
    type?: (string | Array<string> | null);
    properties?: (Record<string, PropertyModel_Input> | null);
    required?: Array<string>;
    description?: (string | null);
    enum?: Array<any>;
    $ref?: any;
    const?: null;
    format?: (string | null);
    default?: null;
    example?: null;
    items?: (ItemsModel_Input | Array<ItemsModel_Input> | null);
};

