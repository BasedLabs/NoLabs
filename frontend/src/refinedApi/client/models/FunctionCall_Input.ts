/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { FunctionCallReturnData_Input } from './FunctionCallReturnData_Input';
import type { FunctionParam } from './FunctionParam';
export type FunctionCall_Input = {
    function_name: string;
    arguments: Array<FunctionParam>;
    data?: FunctionCallReturnData_Input;
};

