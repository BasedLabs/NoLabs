from typing import Any, Dict

from dotenv import load_dotenv

load_dotenv(".env")

from settings import settings
import application
from api_models import SequenceQuery, BlastType
from fastapi import FastAPI
import uvicorn

app = FastAPI(title="Blast")


@app.post("/run_blast")
def run_blast(request: SequenceQuery) -> Dict[str, Any]:
    return application.run_blast(program=BlastType.blastp,
                                 query=request.sequence,
                                 descriptions=request.descriptions,
                                 alignments=request.alignments,
                                 hitlist_size=request.hitlist_size)


uvicorn.run(app, host=settings.fastapi_host, port=settings.fastapi_port)
