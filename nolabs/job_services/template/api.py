from fastapi import FastAPI

app = FastAPI()


@app.post("/run-folding")
def predict_through_api(request: RunEsmFoldPredictionRequest) -> RunEsmFoldPredictionResponse:
    result = run_facebook_api_folding(request)
