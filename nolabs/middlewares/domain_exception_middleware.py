import dataclasses

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse

from nolabs.api_models.problem_details import ProblemDetailsResponse
from nolabs.exceptions import NoLabsException, ErrorCodes


def add_domain_exception_middleware(app: FastAPI):
    @app.middleware("http")
    async def add_process_time_header(request: Request, call_next):
        try:
            response = await call_next(request)
            return response
        except NoLabsException as e:
            if e.error_code == ErrorCodes.experiment_id_not_found:
                return JSONResponse(
                    content={
                        'errors': [e.message],
                        'error_code': e.error_code.value
                    }
                )
            raise e