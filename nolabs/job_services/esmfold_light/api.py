from fastapi import FastAPI

app = FastAPI()


@app.post("/run")
def run(request: PredictFoldingJobRequest) -> PredictFoldingJobResponse:
    result = run_facebook_api_folding(request)
