from fastapi import FastAPI
from p2rank.services import run_p2rank
from p2rank.api_models import RunP2RankPredictionRequest, RunP2RankPredictionResponse
from p2rank.loggers import Log

app = FastAPI()

@app.post("/run-umol")
async def predict(request: RunP2RankPredictionRequest) -> RunP2RankPredictionResponse:
    Log.p2rank_request(request)
    result = run_p2rank(request)
    Log.p2rank_response(result)
    return result

