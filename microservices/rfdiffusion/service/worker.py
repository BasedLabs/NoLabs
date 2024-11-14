from dotenv import load_dotenv

load_dotenv(".env")

import uuid
from typing import Any, Dict

from api_models import RunRfdiffusionRequest
from celery import Celery

from log import logger
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
    broker_transport_options={"heartbeat": 10},
)


@app.task(time_limit=10, name="design")
def design(param: Dict[str, Any]) -> Dict[str, Any]:
    import application

    result = application.design(request=RunRfdiffusionRequest(**param))
    return result.model_dump()


logger.info("Starting celery")
app.worker_main(
    [
        "worker",
        f"--concurrency={settings.celery_worker_concurrency}",
        "-E",
        "-n",
        f"esmfold-{str(uuid.uuid4())}",
    ]
)