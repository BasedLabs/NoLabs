import asyncio

import socketio
from dotenv import load_dotenv

from infrastructure.logging import logger

load_dotenv('infrastructure/.env')

from fastapi import FastAPI, WebSocket
from fastapi.middleware.cors import CORSMiddleware
from starlette.websockets import WebSocketDisconnect

from nolabs.application.middlewares.domain_exception_middleware import add_domain_exception_middleware
from nolabs.application.use_cases.experiments.controller import router as experiment_router
from nolabs.application.use_cases.jobs.controller import router as job_router
from nolabs.application.use_cases.localisation.controller import router as localisation_router
from nolabs.application.use_cases.folding.controller import router as folding_router
from nolabs.application.use_cases.gene_ontology.controller import router as gene_ontology_router
from nolabs.application.use_cases.solubility.controller import router as solubility_router
from nolabs.application.use_cases.conformations.controller import router as conformations_controller
from nolabs.application.use_cases.protein_design.controller import router as protein_design_controller
from nolabs.application.use_cases.binding_pockets.controller import router as binding_pockets_controller
from nolabs.application.use_cases.biobuddy.controller import router as biobuddy_controller
from nolabs.application.use_cases.msa_generation.controller import router as msa_generation_controller
from nolabs.application.use_cases.small_molecules_design.controller import \
    router as small_molecules_design_router
from nolabs.application.event_handlers.di import EventHandlersDependencies
from nolabs.application.use_cases.diffdock.controller import router as diffdock_router
from nolabs.application.use_cases.proteins.controller import router as proteins_router
from nolabs.application.use_cases.ligands.controller import router as ligand_router
from nolabs.infrastructure.logging import setup_logger
from nolabs.application.use_cases.workflow.controller import router as workflow_router
from nolabs.application.use_cases.blast.controller import router as blast_router
from nolabs.infrastructure.mongo_connector import mongo_connect
from nolabs.infrastructure.settings import settings
from nolabs.infrastructure.websocket_queue import websockets_queue

app = FastAPI(
    title='NoLabs',
    version='2.1.7'
)

origins = [
    '*'
]


@app.on_event("startup")
async def startup_event():
    mongo_connect(settings.connection_string)
    EventHandlersDependencies.inject()


sio = socketio.AsyncServer(cors_allowed_origins='*',
                           async_mode='asgi',
                           client_manager=socketio.AsyncRedisManager(settings.socketio_broker))
socket_app = socketio.ASGIApp(sio)

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
app.mount("/socketio", socket_app)

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)

logger.info('Go to /docs to see Swagger')
