from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

import nolabs.controllers.conformations as conformations
import nolabs.infrastructure.environment

app = FastAPI(
    title='NoLabs',
    root_path='/api/v1',
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
app.include_router(conformations.router)