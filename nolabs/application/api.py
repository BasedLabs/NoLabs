from multiprocessing import Process

from dotenv import load_dotenv

load_dotenv("infrastructure/.env")

from nolabs.application.biobuddy.controller import router as biobuddy_controller
from nolabs.application.experiments.controller import router as experiment_router
from nolabs.application.folding.controller import router as folding_router
from nolabs.application.jobs.controller import router as job_router
from nolabs.application.middlewares.domain_exception_middleware import (
    add_domain_exception_middleware,
)
from nolabs.application.proteins.controller import router as proteins_router
from nolabs.workflow.application.controller import router as workflow_router

from nolabs.infrastructure.settings import init_settings
from nolabs.infrastructure.log import init_logging, nolabs_logger as logger
import uvicorn

import socketio
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

settings = init_settings()
init_logging()

app = FastAPI(title="NoLabs")

origins = ["*"]

sio = socketio.AsyncServer(
    cors_allowed_origins="*",
    async_mode="asgi",
    client_manager=socketio.AsyncRedisManager(settings.socketio_broker),
)
socket_app = socketio.ASGIApp(sio)


@sio.event
async def join_room(sid, data):
    experiment_id = data.get("experiment_id")
    if experiment_id:
        await sio.enter_room(sid, experiment_id)
        logger.info(
            "Client joined experiment",
            extra={"client_id": sid, "experiment_id": experiment_id},
        )


app.include_router(experiment_router)
app.include_router(folding_router)
app.include_router(job_router)
app.include_router(proteins_router)
app.include_router(workflow_router)
app.include_router(biobuddy_controller)
add_domain_exception_middleware(app)
app.mount("/", socket_app)

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

logger.info("Go to /docs to see Swagger")

uvicorn.run(app)
