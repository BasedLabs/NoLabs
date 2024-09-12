from fastapi import FastAPI

import application
from api_models import InferenceInput, InferenceOutput

app = FastAPI(title='Esmfold Light')


@app.post("/inference")
def inference(request: InferenceInput) -> InferenceOutput:
    return application.inference(param=request)



