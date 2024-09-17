from fastapi import FastAPI
from sc_gpt.api_models import (EmbedRequest, EmbedResponse,
                               ReferenceMappingRequest,
                               ReferenceMappingResponse)
from sc_gpt.services import classify_cells, get_embeddings

app = FastAPI(title="SC-GPT API")


@app.post("/reference-classify")
def reference_classify(request: ReferenceMappingRequest) -> ReferenceMappingResponse:
    return classify_cells(request)


@app.post("/embed")
def embed(request: EmbedRequest) -> EmbedResponse:
    return get_embeddings(request)
