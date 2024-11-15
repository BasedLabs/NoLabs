from dotenv import load_dotenv

load_dotenv(".env")

from settings import settings

import uvicorn

import application
from api_models import InferenceInput, InferenceOutput
from fastapi import FastAPI

app = FastAPI(title="Esmfold Light")


@app.post("/inference")
def inference(request: InferenceInput) -> InferenceOutput:
    return application.inference(param=request)


uvicorn.run(app, host=settings.fastapi_host, port=settings.fastapi_port)
