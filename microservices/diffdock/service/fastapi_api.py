import application
from api_models import RunDiffDockPredictionRequest, RunDiffDockPredictionResponse
from fastapi import FastAPI

app = FastAPI(title="Diffdock")


@app.post("/inference")
def inference(request: RunDiffDockPredictionRequest) -> RunDiffDockPredictionResponse:
    return application.run_docking(request=request)
