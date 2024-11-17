import uuid
from typing import Any, Dict
from dotenv import load_dotenv

load_dotenv(".env")

from api_models import RunDiffDockPredictionRequest
from celery import Celery
from settings import settings


app = Celery(
    __name__, backend=settings.celery_backend_url, broker=settings.celery_broker_url
)
app.conf.update(
    task_default_queue=settings.celery_worker_queue,
    task_track_started=True,
    task_acks_late=True,
    task_reject_on_worker_lost=True,
    worker_state_db="/opt/celery-state.db",
    accept_content=["application/json", "application/x-python-serialize"],
    broker_transport_options={"heartbeat": 10},
    broker_connection_retry_on_startup=True,
    worker_send_task_events=True
)


@app.task(time_limit=1000, name="inference")
def inference(param: Dict[str, Any]) -> Dict[str, Any]:
    import application

    result = application.run_docking(request=RunDiffDockPredictionRequest(**param))
    return result.model_dump()

app.worker_main(
        [
            "worker",
            f"--concurrency={settings.celery_worker_concurrency}",
            "-E",
            "-n",
            f"diffdock-{str(uuid.uuid4())}",
            "--loglevel", settings.logging_level
        ]
    )
