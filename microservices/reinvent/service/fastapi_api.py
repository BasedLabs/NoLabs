from api_models import RunReinforcementLearningRequest, RunSamplingRequest
from application import Reinvent
from fastapi import APIRouter, FastAPI, File, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from starlette.responses import FileResponse

app = FastAPI(title="Reinvent4 API")

reinvent_router = APIRouter(prefix="/api/reinvent", tags=["reinvent"])


@reinvent_router.post("/start-sampling")
async def sampling(request: RunSamplingRequest):
    """
    Generate new ligands based on the provided config id.
    """
    reinvent_instance = Reinvent()
    return await reinvent_instance.run_sampling_and_scoring(
        config_id=request.config_id,
        number_of_molecules_to_generate=request.number_of_molecules_to_generate,
    )


@reinvent_router.post("/start-learning")
async def learning(request: RunReinforcementLearningRequest):
    """
    Start model learning.
    """
    reinvent_instance = Reinvent()
    return await reinvent_instance.run_reinforcement_learning(request.config_id)


@reinvent_router.post("/prepare-target")
async def prepare_target(pdb_content: UploadFile = File()) -> FileResponse | bytes:
    """
    Prepare .pdbqt file from pdb file.
    """
    reinvent_instance = Reinvent()
    b = await reinvent_instance.prepare_pdbqt(await pdb_content.read())
    return b[0]


origins = ["*"]

app.include_router(reinvent_router)

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
