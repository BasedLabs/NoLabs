import asyncio
from typing import Any, Dict

from api_models import (
    PreparePdbqtRequest,
    PreparePdbqtResponse,
    RunReinforcementLearningRequest,
    RunSamplingRequest,
)
from celery import Celery
from settings import settings

app = Celery(
    __name__, backend=settings.celery_backend_url, broker=settings.celery_broker_url
)
app.conf.update(
    enable_utc=True,
    task_default_queue=settings.celery_worker_queue,
    task_track_started=True,
    task_acks_late=True,
    task_reject_on_worker_lost=True,
    worker_state_db="./celery-state.db",
    accept_content=["application/json", "application/x-python-serialize"]
)


@app.task(time_limit=3 * 24 * 60 * 60, name="reinvent.run_reinforcement_learning")
def run_reinforcement_learning(inp: Dict[str, Any]):
    import application

    reinvent = application.Reinvent()
    request = RunReinforcementLearningRequest(**inp)
    reinvent.run_reinforcement_learning(config_id=request.config_id)


@app.task(time_limit=3 * 24 * 60 * 60, name="reinvent.run_sampling")
def run_sampling(inp: Dict[str, Any]):
    import application

    reinvent = application.Reinvent()
    request = RunSamplingRequest(**inp)
    reinvent.run_sampling_and_scoring(
        config_id=request.config_id,
        number_of_molecules_to_generate=request.number_of_molecules_to_generate,
    )


@app.task(time_limit=30, name="reinvent.prepare_target")
def prepare_target(inp: Dict[str, Any]):
    import application

    reinvent = application.Reinvent()
    request = PreparePdbqtRequest(**inp)
    result = asyncio.run(reinvent.prepare_pdbqt(pdb=request.pdb))
    return PreparePdbqtResponse(pdbqt=result[0], file_path=result[1]).model_dump()
