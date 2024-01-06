from fastapi import FastAPI
from umol.services import run_umol
from umol.api_models import RunUmolPredictionRequest, RunUmolPredictionResponce
from umol.loggers import Log

app = FastAPI()


@app.post("/run-umol")
async def predict(request: RunUmolPredictionRequest) -> RunUmolPredictionResponce:
    Log.umol_request(request)
    result = run_umol(request)
    Log.umol_response(result)
    return result

