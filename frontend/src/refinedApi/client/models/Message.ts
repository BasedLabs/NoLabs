/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { FunctionCall_Output } from './FunctionCall_Output';
import type { RegularMessage } from './RegularMessage';
export type Message = {
    id: string;
    role: string;
    message: (RegularMessage | Array<FunctionCall_Output>);
    type: string;
};

