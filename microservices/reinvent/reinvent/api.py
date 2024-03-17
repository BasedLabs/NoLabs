import asyncio
from typing import List

from fastapi.middleware.cors import CORSMiddleware

from fastapi import FastAPI, File, APIRouter, UploadFile

from reinvent.api_models import RunFineTuningJobRequest, FineTuningJobResponse, RunInferenceRequest, \
    RunInferenceResponse
from reinvent.services import FineTuning, Inference

app = FastAPI(
    title='Reinvent4 API',
    version='0.1'
)

ft_router = APIRouter(
    prefix='/api/v1/fine-tuning',
    tags=['folding']
)

inference_router = APIRouter(
    prefix='/api/v1/inference',
    tags=['inference']
)

FineTuning().run_refresh_progress_task()

@ft_router.post("/run")
async def run_fine_tuning(request: RunFineTuningJobRequest, pdb_content: UploadFile = File()) -> FineTuningJobResponse:
    ft = FineTuning()
    return await ft.run(pdb_content, request)


@ft_router.get("/{job_id}")
async def fine_tuning_get_job(job_id: str) -> FineTuningJobResponse | None:
    ft = FineTuning()
    return ft.get_job(job_id)


@ft_router.get("/all")
async def fine_tuning_get_all_jobs() -> List[FineTuningJobResponse]:
    ft = FineTuning()
    return ft.all_jobs()


@inference_router.post("/run")
async def run_inference(request: RunInferenceRequest) -> RunInferenceResponse:
    inference = Inference()
    return inference.run(request)


origins = [
    '*'
]

app.include_router(ft_router)

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)

print('Go to /api/v1/docs to see Swagger')
