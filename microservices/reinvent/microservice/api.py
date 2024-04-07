from typing import List, Optional

from fastapi.middleware.cors import CORSMiddleware

from fastapi import FastAPI, File, UploadFile
from starlette.responses import FileResponse

from microservice.api_models import JobResponse, ParamsResponse, LogsResponse, SmilesResponse, \
    ParamsRequest
from microservice.services import Inference, RunType

app = FastAPI(
    title='Reinvent4 API',
    version='0.1'
)


@app.post("/jobs/{job_id}/sampling")
async def sampling(
        job_id: str):
    inference = Inference()
    return await inference.run(job_id, RunType.SAMPLING_SCORING)


@app.post("/jobs/{job_id}/learning")
async def learning(
        job_id: str):
    inference = Inference()
    return await inference.run(job_id, RunType.RL)


@app.get('/jobs/{job_id}/params')
async def params(job_id: str) -> Optional[ParamsResponse]:
    inference = Inference()
    return inference.get_params(job_id)


@app.post('/jobs/{job_id}/params')
async def save_params(job_id: str,
                      name: str,
                      center_x: float,
                      center_y: float, center_z: float, size_x: float,
                      size_y: float, size_z: float,
                      epochs: int = 50,
                      batch_size: int = 128,
                      minscore: float = 0.4,
                      pdb_file: UploadFile = File()):
    inference = Inference()
    await inference.save_params(pdb_file=pdb_file,
                   request=ParamsRequest(job_id=job_id,name=name, center_x=center_x, center_y=center_y, center_z=center_z,
                                         size_x=size_x,
                                         size_y=size_y, size_z=size_z, batch_size=batch_size, minscore=minscore,
                                         epochs=epochs))


@app.post("/jobs/{job_id}/stop")
async def stop(
        job_id: str):
    inference = Inference()
    return await inference.stop(job_id)


@app.delete("/jobs/{job_id}")
async def delete(
        job_id: str):
    inference = Inference()
    return await inference.delete(job_id)


@app.post('/prepare-binder')
async def prepare_binder(pdb_content: UploadFile = File()) -> FileResponse:
    inference = Inference()
    return await inference.prepare_pdbqt(pdb_content)


@app.get("/jobs/{job_id}")
async def get_job(job_id: str) -> Optional[JobResponse]:
    inference = Inference()
    return inference.get_job(job_id)


@app.get("/jobs")
async def get_all_jobs() -> List[JobResponse]:
    inference = Inference()
    return inference.all_jobs()


@app.get('/jobs/{job_id}/logs')
async def logs(job_id: str) -> Optional[LogsResponse]:
    inference = Inference()
    return inference.get_logs(job_id)


@app.get('/jobs/{job_id}/smiles')
async def smiles(job_id: str) -> SmilesResponse:
    inference = Inference()
    return inference.get_smiles(job_id)


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
