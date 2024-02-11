import {ErrorResponse} from "src/api/errorTypes";

export const obtainErrorResponse = (response: any): ErrorResponse | null => {
    if(response && 'errors' in response && 'error_code' in response){
        return response as ErrorResponse;
    }
    return null;
}
