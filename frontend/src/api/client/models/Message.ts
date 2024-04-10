/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { FunctionCall } from './FunctionCall';
import type { RegularMessage } from './RegularMessage';
export type Message = {
    role: string;
    message: (RegularMessage | Array<FunctionCall>);
    type: string;
};

