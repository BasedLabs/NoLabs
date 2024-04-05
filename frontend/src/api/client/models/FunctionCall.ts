/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { FunctionCallReturnData } from './FunctionCallReturnData';
import type { FunctionParam } from './FunctionParam';
export type FunctionCall = {
    function_name: string;
    parameters: Array<FunctionParam>;
    data?: FunctionCallReturnData;
};

