from dotenv import load_dotenv

load_dotenv(".env")

from log import logger
from settings import settings
import application
from api_models import RunRfdiffusionRequest, RunRfdiffusionResponse
from fastapi import FastAPI

import uvicorn

app = FastAPI(title="Esmfold")


@app.post("/design")
def inference(request: RunRfdiffusionRequest) -> RunRfdiffusionResponse:
    return application.design(request=request)


logger.info("Starting fastapi")
uvicorn.run(app, host=settings.fastapi_host, port=settings.fastapi_port)
