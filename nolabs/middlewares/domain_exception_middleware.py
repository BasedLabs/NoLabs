from fastapi import FastAPI, Request, HTTPException
from fastapi.responses import JSONResponse

from nolabs.exceptions import NoLabsException, ErrorCodes


def add_domain_exception_middleware(app: FastAPI):
    @app.middleware("http")
    async def add_process_time_header(request: Request, call_next):
        try:
            response = await call_next(request)
            return response
        except NoLabsException as e:
            return JSONResponse(content={
                'errors': e.messages,
                'error_code': e.error_code.value
            }, headers={'Content-Type': 'application/problem+json'}, status_code=200)
        except Exception as e:
            print(str(e))
            return JSONResponse(content={
                'errors': ['Unknown server error'],
                'error_code': ErrorCodes.unknown_exception.value,
            }, headers={'Content-Type': 'application/problem+json'},
                status_code=200)
