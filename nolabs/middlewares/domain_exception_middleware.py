from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.infrastructure.nolabs_logging import logger


def add_domain_exception_middleware(app: FastAPI):
    @app.middleware("http")
    async def add_process_time_header(request: Request, call_next):
        try:
            response = await call_next(request)
            return response
        except NoLabsException as e:
            return JSONResponse(content={
                'errors': e.messages,
                'error_code': e.error_code
            }, headers={'Content-Type': 'application/problem+json'}, status_code=200)
        except Exception as e:
            logger.exception('Exception occured in application')
            return JSONResponse(content={
                'errors': ['Unknown server error'],
                'error_code': ErrorCodes.unknown_exception.value,
            }, headers={'Content-Type': 'application/problem+json'},
                status_code=200)
