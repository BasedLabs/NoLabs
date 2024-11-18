from dotenv import load_dotenv

load_dotenv(".env")

from settings import settings
import application
from api_models import RunRfdiffusionRequest, RunRfdiffusionResponse
from fastapi import FastAPI

import uvicorn

app = FastAPI(title="Rfdiffusion")


@app.post("/design")
def inference(request: RunRfdiffusionRequest) -> RunRfdiffusionResponse:
    return application.design(request=request)


uvicorn.run(app, host=settings.fastapi_host, port=settings.fastapi_port)
