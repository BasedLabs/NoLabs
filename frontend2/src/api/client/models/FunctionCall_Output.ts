/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { FunctionCallReturnData_Output } from './FunctionCallReturnData_Output';
import type { FunctionParam } from './FunctionParam';
export type FunctionCall_Output = {
    function_name: string;
    arguments: Array<FunctionParam>;
    data?: FunctionCallReturnData_Output;
};

