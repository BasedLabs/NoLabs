from typing import Any, Dict

from celery import Celery

from api_models import RunDiffDockPredictionRequest
from settings import settings

print('hey')

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


@app.task(time_limit=600, name="diffdock.inference")
def inference(param: Dict[str, Any]) -> Dict[str, Any]:
    import application

    result = application.run_docking(request=RunDiffDockPredictionRequest(**param))
    return result.model_dump()
