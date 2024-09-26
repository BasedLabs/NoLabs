from typing import Any, Dict

from api_models import InferenceInput
from celery import Celery
from settings import settings

app = Celery(
    __name__, backend=settings.celery_backend_url, broker=settings.celery_broker_url
)
app.conf.enable_utc = settings.celery_enable_utc
app.conf.task_default_queue = settings.celery_worker_queue
app.conf.accept_content = ["application/json", "application/x-python-serialize"]


@app.task(time_limit=10, name="esmfold-light-service.inference")
def inference(param: Dict[str, Any]) -> Dict[str, Any]:
    import application

    result = application.inference(param=InferenceInput(**param))
    return result.dict()
