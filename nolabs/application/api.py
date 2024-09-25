from contextlib import asynccontextmanager

import uvicorn
from dotenv import load_dotenv

load_dotenv("infrastructure/.env")

import socketio
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from nolabs.application.biobuddy.controller import \
    router as biobuddy_controller
from nolabs.application.diffdock.controller import router as diffdock_router
from nolabs.application.event_handlers.di import EventHandlersDependencies
from nolabs.application.experiments.controller import \
    router as experiment_router
from nolabs.application.folding.controller import router as folding_router
from nolabs.application.ligands.controller import router as ligand_router
from nolabs.application.middlewares.domain_exception_middleware import \
    add_domain_exception_middleware
from nolabs.application.proteins.controller import router as proteins_router
from nolabs.application.small_molecules_design.controller import \
    router as small_molecules_design_router
from nolabs.application.use_cases.binding_pockets.controller import \
    router as binding_pockets_controller
from nolabs.application.use_cases.blast.controller import \
    router as blast_router
from nolabs.application.use_cases.conformations.controller import \
    router as conformations_controller
from nolabs.application.use_cases.gene_ontology.controller import \
    router as gene_ontology_router
from nolabs.application.use_cases.jobs.controller import router as job_router
from nolabs.application.use_cases.localisation.controller import \
    router as localisation_router
from nolabs.application.use_cases.msa_generation.controller import \
    router as msa_generation_controller
from nolabs.application.use_cases.protein_design.controller import \
    router as protein_design_controller
from nolabs.application.use_cases.solubility.controller import \
    router as solubility_router
from nolabs.application.workflow.api.controller import \
    router as workflow_router
from nolabs.infrastructure.log import logger
from nolabs.infrastructure.mongo_connector import mongo_connect
from nolabs.infrastructure.settings import settings


@asynccontextmanager
async def lifespan(app: FastAPI):
    mongo_connect(settings.connection_string)
    EventHandlersDependencies.inject()
    yield


app = FastAPI(title="NoLabs", version="2.1.7", lifespan=lifespan)

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


app.include_router(localisation_router)
app.include_router(experiment_router)
app.include_router(folding_router)
app.include_router(job_router)
app.include_router(gene_ontology_router)
app.include_router(solubility_router)
app.include_router(conformations_controller)
app.include_router(protein_design_controller)
app.include_router(binding_pockets_controller)
app.include_router(msa_generation_controller)
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
