/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { FunctionCallReturnData } from './FunctionCallReturnData';
import type { FunctionParam } from './FunctionParam';
export type FunctionCall = {
    function_name: string;
    arguments: Array<FunctionParam>;
    data?: FunctionCallReturnData;
};

