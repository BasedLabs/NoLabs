from dotenv import load_dotenv

load_dotenv(".env")

import uvicorn
from fastapi import FastAPI

from settings import settings
from application import run_protmpnn
from api_models import RunProtMPNNPredictionRequest, RunProtMPNNPredictionResponse

app = FastAPI(
    title="ProtMPNN"
)

@app.post("/run-design")
def predict(request: RunProtMPNNPredictionRequest) -> RunProtMPNNPredictionResponse:
    return run_protmpnn(request)

uvicorn.run(app, host=settings.fastapi_host, port=settings.fastapi_port)