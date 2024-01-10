import uvicorn
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from nolabs.controllers.conformations.conformations import router as conformations_router
from nolabs.controllers.solubility.solubility import router as solubility_router
import nolabs.infrastructure.environment

pfx = '/api/v1'

app = FastAPI(
    title='NoLabs',
    version='1'
)

origins = [
    '*'
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)
app.include_router(conformations_router)
app.include_router(solubility_router)

print('Go to /api/v1/docs to see Swagger')