from fastapi import FastAPI
from gene_ontology.api_models import (UUID, IsJobRunningResponse,
                                      RunGeneOntologyPredictionRequest,
                                      RunGeneOntologyPredictionResponse)
from gene_ontology.services import run_gene_ontology_prediction
from gene_ontology.shared import JobStatusEnum, memory_manager

app = FastAPI(title="Gene ontology api")

from gene_ontology.loggers import logger

logger.starting_api()


@app.post("/run")
async def run_go_prediction(
    request: RunGeneOntologyPredictionRequest,
) -> RunGeneOntologyPredictionResponse:
    logger.run_go_prediction_request(request)
    result = run_gene_ontology_prediction(request)
    logger.run_go_prediction_response(result)
    return result


@app.get("/job/{job_id}/is-running")
def is_job_running(job_id: UUID) -> IsJobRunningResponse:
    return IsJobRunningResponse(
        is_running=memory_manager.get_job_status(str(job_id)) == JobStatusEnum.running
    )
