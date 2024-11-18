from dotenv import load_dotenv
from settings import settings

load_dotenv(".env")

import uvicorn
import application
from api_models import RunDiffDockPredictionRequest, RunDiffDockPredictionResponse
from fastapi import FastAPI

app = FastAPI(title="Diffdock")


@app.post("/inference")
def inference(request: RunDiffDockPredictionRequest) -> RunDiffDockPredictionResponse:
    return application.run_docking(request=request)



uvicorn.run(app, host=settings.fastapi_host, port=settings.fastapi_port)
