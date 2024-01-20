from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from nolabs.controllers.conformations.conformations import router as conformations_router
from nolabs.controllers.solubility.solubility import router as solubility_router
from nolabs.controllers.localisation.localisation import router as localisation_router
from nolabs.controllers.gene_ontology.gene_ontology import router as gene_ontology_router
from nolabs.controllers.drug_discovery.drug_discovery import router as drug_discovery_router
from nolabs.controllers.protein_design.protein_design import router as protein_design_router
from nolabs.middlewares.domain_exception_middleware import add_domain_exception_middleware
import nolabs.infrastructure.environment

pfx = '/api/v1'

app = FastAPI(
    title='NoLabs',
    version='1'
)

origins = [
    '*'
]

add_domain_exception_middleware(app)
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)

app.include_router(conformations_router)
app.include_router(solubility_router)
app.include_router(localisation_router)
app.include_router(gene_ontology_router)
app.include_router(drug_discovery_router)
app.include_router(protein_design_router)

print('Go to /api/v1/docs to see Swagger')