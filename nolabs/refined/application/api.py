import sys

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from nolabs.middlewares.domain_exception_middleware import add_domain_exception_middleware
from nolabs.refined.application.controllers.amino_acid.localisation.controller import router as localisation_router
from nolabs.refined.application.controllers.experiments.controller import router as experiment_router
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
add_domain_exception_middleware(app)

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)

print('Go to /api/v1/docs to see Swagger')
