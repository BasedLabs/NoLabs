from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse

from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.infrastructure.log import logger


def add_domain_exception_middleware(app: FastAPI):
    @app.exception_handler(NoLabsException)
    async def handle_nolabs_exception(request: Request, exc: NoLabsException):
        logger.exception(
            exc.message,
            exc_info=exc,
            extra={"httpRequest": {"requestUrl": request.url}},
        )
        return JSONResponse(
            content={
                "message": exc.message,
                "error_code": exc.error_code,
                "data": exc.data,
            },
            headers={"Content-Type": "application/problem+json"},
            status_code=500,
        )

    @app.exception_handler(Exception)
    async def handle_exception(request: Request, exc: Exception):
        logger.exception(
            str(exc),
            exc_info=exc,
            extra={"httpRequest": {"requestUrl": request.url}},
        )
        return JSONResponse(
            content={
                "errors": ["Unknown server error"],
                "error_code": ErrorCodes.unknown_exception.value.code,
            },
            headers={"Content-Type": "application/problem+json"},
            status_code=500,
        )
