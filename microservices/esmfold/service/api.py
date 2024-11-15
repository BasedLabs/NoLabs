from dotenv import load_dotenv

load_dotenv(".env")

from settings import settings
import application
from api_models import InferenceInput, InferenceOutput
from fastapi import FastAPI


import uvicorn

app = FastAPI(title="Esmfold")


@app.post("/inference")
def inference(request: InferenceInput) -> InferenceOutput:
    return application.inference(param=request)

uvicorn.run(app, host=settings.fastapi_host, port=settings.fastapi_port)
