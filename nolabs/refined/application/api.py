import logging

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from nolabs.refined.application.middlewares.domain_exception_middleware import add_domain_exception_middleware
from nolabs.refined.application.use_cases.localisation.controller import router as localisation_router
from nolabs.refined.application.use_cases.experiments.controller import router as experiment_router
from nolabs.refined.application.use_cases.jobs.controller import router as job_router
from nolabs.refined.application.use_cases.folding.controller import router as folding_router
from nolabs.refined.application.use_cases.gene_ontology.controller import router as gene_ontology_router
from nolabs.refined.application.use_cases.solubility.controller import router as solubility_router
from nolabs.refined.application.use_cases.conformations.controller import router as conformations_controller
from nolabs.refined.application.use_cases.protein_design.controller import router as protein_design_controller
from nolabs.refined.application.use_cases.binding_pockets.controller import router as binding_pockets_controller
from nolabs.refined.application.use_cases.biobuddy.controller import router as biobuddy_controller
from nolabs.refined.application.use_cases.msa_generation.controller import router as msa_generation_controller
from nolabs.refined.application.use_cases.umol.controller import router as umol_controller
from nolabs.refined.application.use_cases.small_molecules_design.controller import \
    router as small_molecules_design_router
from nolabs.refined.application.event_handlers.di import EventHandlersDependencies
from nolabs.refined.application.use_cases.diffdock.controller import router as diffdock_router
from nolabs.refined.application.use_cases.proteins.controller import router as proteins_router
from nolabs.refined.application.use_cases.ligands.controller import router as ligand_router
from nolabs.refined.infrastructure.logging import setup_logger
from nolabs.workflow.application.controller import router as workflow_router
from nolabs.refined.infrastructure.mongo_connector import mongo_connect
from nolabs.refined.infrastructure.settings import Settings

app = FastAPI(
    title='NoLabs',
    version='2.0.0'
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
app.include_router(solubility_router)
app.include_router(conformations_controller)
app.include_router(protein_design_controller)
app.include_router(binding_pockets_controller)
app.include_router(msa_generation_controller)
app.include_router(small_molecules_design_router)
app.include_router(diffdock_router)
app.include_router(umol_controller)
app.include_router(proteins_router)
app.include_router(ligand_router)
app.include_router(workflow_router)
app.include_router(biobuddy_controller)
add_domain_exception_middleware(app)

setup_logger()


app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)

print('Go to /api/v1/docs to see Swagger')
