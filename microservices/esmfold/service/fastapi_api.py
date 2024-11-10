import application
from api_models import InferenceInput, InferenceOutput
from fastapi import FastAPI

app = FastAPI(title="Esmfold")


@app.post("/inference")
def inference(request: InferenceInput) -> InferenceOutput:
    return application.inference(param=request)
