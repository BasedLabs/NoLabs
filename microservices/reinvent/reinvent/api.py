from typing import List, Optional

from fastapi.middleware.cors import CORSMiddleware

from fastapi import FastAPI, File, UploadFile, APIRouter
from starlette.responses import FileResponse

from reinvent.api_models import ConfigurationResponse, ParamsResponse, LogsResponse, SmilesResponse, \
    ParamsRequest
from reinvent.services import Reinvent, RunType

from reinvent.api_models import SamplingSizeRequest

app = FastAPI(
    title='Reinvent4 API',
    version='0.0.1'
)

reinvent_router = APIRouter(prefix='/api/reinvent', tags=['reinvent'])
preparation_router = APIRouter(prefix='/api/preparation', tags=['preparation'])


@reinvent_router.post("/{config_id}/start-sampling")
async def sampling(
        config_id: str,
        request: SamplingSizeRequest):
    """
    Generate new ligands based on the provided config id.
    """
    reinvent_instance = Reinvent()
    return await reinvent_instance.run(config_id, RunType.SAMPLING_SCORING, sampling_size=request.number_of_molecules_to_design)


@reinvent_router.post("/{config_id}/start-learning")
async def learning(
        config_id: str):
    """
    Start model learning.
    """
    reinvent_instance = Reinvent()
    return await reinvent_instance.run(config_id, RunType.RL)


@reinvent_router.get('/{config_id}/params')
async def params(config_id: str) -> Optional[ParamsResponse]:
    """
    Get learning parameters of configuration.
    """
    reinvent_instance = Reinvent()
    return reinvent_instance.get_params(config_id)


@reinvent_router.post('/{config_id}/params')
async def save_params(config_id: str,
                      name: str,
                      center_x: float,
                      center_y: float, center_z: float, size_x: float,
                      size_y: float, size_z: float,
                      epochs: int = 50,
                      batch_size: int = 128,
                      minscore: float = 0.4,
                      pdb_file: UploadFile = File()):
    """
    Save parameters for reinvent reinforcement learning configuration.
    """
    reinvent_instance = Reinvent()
    await reinvent_instance.save_params(pdb_file=pdb_file,
                   request=ParamsRequest(config_id=config_id,name=name, center_x=center_x, center_y=center_y, center_z=center_z,
                                         size_x=size_x,
                                         size_y=size_y, size_z=size_z, batch_size=batch_size, minscore=minscore,
                                         epochs=epochs))


@reinvent_router.post("/{config_id}/jobs/stop")
async def stop(
        config_id: str):
    """
    Stop current job.
    """
    reinvent_instance = Reinvent()
    return await reinvent_instance.stop(config_id)


@reinvent_router.delete("/{config_id}")
async def delete(
        config_id: str):
    """
    Delete configuration.
    """
    reinvent_instance = Reinvent()
    return await reinvent_instance.delete(config_id)


@preparation_router.post('/prepare-target')
async def prepare_target(pdb_content: UploadFile = File()) -> FileResponse:
    """
    Prepare .pdbqt file from pdb file.
    """
    reinvent_instance = Reinvent()
    return await reinvent_instance.prepare_pdbqt(pdb_content)


@reinvent_router.get("/reinvent/{config_id}")
async def get_config(config_id: str) -> Optional[ConfigurationResponse]:
    """
    Get configuration.
    """
    reinvent_instance = Reinvent()
    return reinvent_instance.get_config(config_id)


@reinvent_router.get("/")
async def get_all_configs() -> List[ConfigurationResponse]:
    """
    Get all configurations available.
    """
    reinvent_instance = Reinvent()
    return reinvent_instance.all_configs()


@reinvent_router.get('/{config_id}/logs')
async def logs(config_id: str) -> Optional[LogsResponse]:
    """
    Get logs of all runs jobs of given configuration.
    """
    reinvent_instance = Reinvent()
    return reinvent_instance.get_logs(config_id)


@reinvent_router.get('/{config_id}/smiles')
async def smiles(config_id: str) -> SmilesResponse:
    """
    Get generated smiles after sampling or RL.
    """
    reinvent_instance = Reinvent()
    return reinvent_instance.get_smiles(config_id)


origins = [
    '*'
]

app.include_router(reinvent_router)
app.include_router(preparation_router)

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)

print('/docs - Swagger')
