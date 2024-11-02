from typing import Any, Dict

from api_models import InferenceInput
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
    worker_state_db="/opt/celery-state.db",
    accept_content=["application/json", "application/x-python-serialize"],
    broker_transport_options={'heartbeat': 10}
)


@app.task(time_limit=10, name="inference")
def inference(param: Dict[str, Any]) -> Dict[str, Any]:
    import application

    result = application.inference(param=InferenceInput(**param))
    return result.dict()
