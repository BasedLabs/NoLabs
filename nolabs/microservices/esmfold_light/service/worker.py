import logging
from typing import Any, Dict

from settings import settings

from celery import Celery

logger = logging.getLogger(__name__)

logger.info('Starting celery service')

from api_models import InferenceInput

app = Celery(__name__, backend=settings.celery_backend_url, broker=settings.celery_broker_url)
app.conf.enable_utc = settings.celery_enable_utc
app.conf.accept_content = ['application/json', 'application/x-python-serialize']


@app.task(time_limit=10, name='esmfold-light-service.inference')
def inference(param: Dict[str, Any]) -> Dict[str, Any]:
    import application
    result = application.inference(param=InferenceInput(**param))
    return result.dict()
