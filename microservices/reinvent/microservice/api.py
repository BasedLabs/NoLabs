from typing import List

from fastapi.middleware.cors import CORSMiddleware

from fastapi import FastAPI, File, APIRouter, UploadFile
from starlette.responses import FileResponse

from microservice.api_models import RunFineTuningJobRequest, FineTuningJobResponse, RunInferenceRequest, \
    RunInferenceResponse
from microservice.services import FineTuning

app = FastAPI(
    title='Reinvent4 API',
    version='0.1'
)

ft_router = APIRouter(
    prefix='/api/v1/fine-tuning',
    tags=['fine-tuning']
)

inference_router = APIRouter(
    prefix='/api/v1/inference',
    tags=['inference']
)


@ft_router.post("/run")
async def run_fine_tuning(center_x: float,
                          center_y: float, center_z: float, size_x: float,
                          size_y: float, size_z: float,
                          epochs: int = 50,
                          batch_size: int = 128,
                          minscore: float = 0.4,
                          pdb_content: UploadFile = File()) -> FineTuningJobResponse:
    ft = FineTuning()
    return await ft.run(pdb_content,
                        RunFineTuningJobRequest(center_x=center_x, center_y=center_y, center_z=center_z, size_x=size_x,
                                                size_y=size_y, size_z=size_z, batch_size=batch_size, minscore=minscore, epochs=epochs))


@ft_router.post('/prepare-binder')
async def prepare_binder(pdb_content: UploadFile = File()) -> FileResponse:
    ft = FineTuning()
    return await ft.prepare_pdbqt(pdb_content)


@ft_router.get("/{job_id}")
async def fine_tuning_get_job(job_id: str) -> FineTuningJobResponse | None:
    ft = FineTuning()
    return ft.get_job(job_id)


@ft_router.get("/all")
async def fine_tuning_get_all_jobs() -> List[FineTuningJobResponse]:
    ft = FineTuning()
    return ft.all_jobs()


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
