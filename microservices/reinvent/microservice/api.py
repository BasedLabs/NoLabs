from typing import List

from fastapi.middleware.cors import CORSMiddleware

from fastapi import FastAPI, File, APIRouter, UploadFile
from starlette.responses import FileResponse

from microservice.api_models import JobResponse, RunJobRequest, ParamsResponse, LogsResponse, SmilesResponse
from microservice.services import FineTuning

app = FastAPI(
    title='Reinvent4 API',
    version='0.1'
)


@app.post("/jobs/run")
async def run(
        name: str,
        center_x: float,
        center_y: float, center_z: float, size_x: float,
        size_y: float, size_z: float,
        epochs: int = 50,
        batch_size: int = 128,
        minscore: float = 0.4,
        pdb_content: UploadFile = File()) -> JobResponse:
    ft = FineTuning()
    return await ft.run(pdb_content,
                        RunJobRequest(name=name, center_x=center_x, center_y=center_y, center_z=center_z, size_x=size_x,
                                      size_y=size_y, size_z=size_z, batch_size=batch_size, minscore=minscore,
                                      epochs=epochs))


@app.post('/prepare-binder')
async def prepare_binder(pdb_content: UploadFile = File()) -> FileResponse:
    ft = FineTuning()
    return await ft.prepare_pdbqt(pdb_content)


@app.get("/jobs/{job_id}")
async def get_job(job_id: str) -> JobResponse | None:
    ft = FineTuning()
    return ft.get_job(job_id)


@app.get("/jobs")
async def fine_tuning_get_all_jobs() -> List[JobResponse]:
    ft = FineTuning()
    return ft.all_jobs()

@app.get("/jobs/{job_id}/params")
async def params(job_id: str) -> ParamsResponse | None:
    ft = FineTuning()
    return ft.get_params(job_id)

@app.get('/jobs/{job_id}/logs')
async def logs(job_id: str) -> LogsResponse | None:
    ft = FineTuning()
    return ft.get_logs(job_id)

@app.get('/jobs/{job_id}/smiles')
async def smiles(job_id: str) -> SmilesResponse:
    ft = FineTuning()
    return ft.get_smiles(job_id)


origins = [
    '*'
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)

print('Go to /api/v1/docs to see Swagger')
