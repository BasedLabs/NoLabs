from fastapi import FastAPI
from p2rank.services import run_p2rank
from p2rank.api_models import RunP2RankPredictionRequest, RunP2RankPredictionResponse

app = FastAPI()

from p2rank.loggers import logger

@app.post("/run-p2rank")
async def predict(request: RunP2RankPredictionRequest) -> RunP2RankPredictionResponse:
    logger.p2rank_request(request)
    result = run_p2rank(request)
    logger.p2rank_response(result)
    return result

