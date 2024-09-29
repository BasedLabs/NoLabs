from contextlib import asynccontextmanager
from multiprocessing import Process

from dotenv import load_dotenv


load_dotenv("infrastructure/.env")

from nolabs.infrastructure.settings import settings


def serve_fastapi():
    import uvicorn

    import socketio
    from fastapi import FastAPI
    from fastapi.middleware.cors import CORSMiddleware

    from nolabs.application import initialize
    from nolabs.application.biobuddy.controller import router as biobuddy_controller
    from nolabs.application.blast.controller import router as blast_router
    from nolabs.application.diffdock.controller import router as diffdock_router
    from nolabs.application.experiments.controller import router as experiment_router
    from nolabs.application.folding.controller import router as folding_router
    from nolabs.application.jobs.controller import router as job_router
    from nolabs.application.ligands.controller import router as ligand_router
    from nolabs.application.middlewares.domain_exception_middleware import (
        add_domain_exception_middleware,
    )
    from nolabs.application.proteins.controller import router as proteins_router
    from nolabs.application.small_molecules_design.controller import (
        router as small_molecules_design_router,
    )
    from workflow.controller import router as workflow_router
    from nolabs.infrastructure.log import logger
    from nolabs.infrastructure.settings import settings

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        initialize()
        yield

    app = FastAPI(title="NoLabs", lifespan=lifespan)

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
    app.include_router(small_molecules_design_router)
    app.include_router(diffdock_router)
    app.include_router(proteins_router)
    app.include_router(ligand_router)
    app.include_router(workflow_router)
    app.include_router(biobuddy_controller)
    app.include_router(blast_router)
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


def serve_workflow():
    from workflow.flows import serve_workflow

    serve_workflow()


def serve_united():
    Process(target=serve_workflow).start()
    serve_fastapi()


if settings.mode == "fastapi":
    serve_fastapi()
if settings.mode == "workflow":
    serve_workflow()
if settings.mode == "united":
    serve_united()
