from fastapi import FastAPI
app = FastAPI()

from msa_light.api_models import RunMsaPredictionRequest, RunMsaPredictionResponse
from msa_light.loggers import Log
from msa_light.services import predict_msa_service

@app.post("/predict-msa")
def predict_msa(request: RunMsaPredictionRequest) -> RunMsaPredictionResponse:
    Log.msa_request(request)
    result = predict_msa_service(request)
    Log.msa_response(result)
    return result

