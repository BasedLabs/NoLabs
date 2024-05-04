import sys

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from nolabs.middlewares.domain_exception_middleware import add_domain_exception_middleware
from nolabs.refined.application.amino_acid.localisation.controller import router as localisation_router
from nolabs.refined.application.experiments.controller import router as experiment_router
from nolabs.refined.application.jobs.controller import router as job_router
from nolabs.refined.application.amino_acid.folding.controller import router as folding_router
from nolabs.refined.application.amino_acid.gene_ontology.controller import router as gene_ontology_router
from nolabs.refined.application.conformations.controller import router as conformations_controller
from nolabs.refined.application.protein_design.controller import router as protein_design_controller
from nolabs.refined.application.event_handlers.di import EventHandlersDependencies
from nolabs.refined.infrastructure.mongo_connector import mongo_connect
from nolabs.refined.infrastructure.settings import Settings

app = FastAPI(
    title='NoLabs',
    version='1.1.0'
)

origins = [
    '*'
]


@app.on_event("startup")
async def startup_event():
    settings = Settings.load()
    mongo_connect(settings.connection_string)
    EventHandlersDependencies.inject()


app.include_router(localisation_router)
app.include_router(experiment_router)
app.include_router(folding_router)
app.include_router(job_router)
app.include_router(gene_ontology_router)
app.include_router(conformations_controller)
app.include_router(protein_design_controller)
add_domain_exception_middleware(app)

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)

print('Go to /api/v1/docs to see Swagger')