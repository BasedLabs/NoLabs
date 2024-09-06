from celery import Celery
from pydantic import BaseModel

from nolabs.infrastructure.settings import settings

celery_app = Celery(__name__, broker=settings.celery_broker, backend=settings.celery_backend)
celery_app.conf.enable_utc = True


def send_selery_task(name: str, payload: BaseModel):
    return celery_app.send_task(
        name=name,
        args=[payload.dict()])
