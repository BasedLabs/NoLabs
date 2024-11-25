__all__ = [
    "router",
]

import uuid
from typing import List
from uuid import UUID

from fastapi import APIRouter, Request

from nolabs.application.adaptyv_bio.protein_affinity_characterization.api_models import TargetResponse, \
    EstimatesResponse, JobResponse, SetupJobRequest
from nolabs.application.adaptyv_bio.protein_affinity_characterization.use_cases import StartJobFeature, GetJobFeature, \
    SetupJobFeature, ListTargetsFeature, GetEstimatesFeature

router = APIRouter(prefix="/api/v1/adaptyv-bio/protein-affinity", tags=["Protein affinity characterization"])


@router.post("/jobs/run/{job_id}", summary="Run job")
async def start_job(job_id: UUID):
    await StartJobFeature().handle(job_id)


@router.get("/jobs/{job_id}", summary="Get job")
async def get_job(job_id: UUID) -> JobResponse:
    return await GetJobFeature().handle(job_id=job_id)


@router.post("/jobs", summary="Setup job")
async def setup_job(request: SetupJobRequest, http_request: Request) -> JobResponse:
    return await SetupJobFeature().handle(request=request, session_url=http_request.url.path)

@router.get("/list-targets/{search_query}", summary="List targets for the experiment")
async def list_targets(search_query: str) -> List[TargetResponse]:
    return await ListTargetsFeature().handle(search_query=search_query)

@router.get('/jobs/{job_id}/estimates', summary="Get experiment estimates")
async def get_estimates(job_id: uuid.UUID) -> EstimatesResponse:
    return await GetEstimatesFeature().handle(job_id=job_id)