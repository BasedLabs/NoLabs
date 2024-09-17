import asyncio
from typing import Any

from celery import Celery
from celery.result import AsyncResult
from pydantic import BaseModel

from nolabs.infrastructure.settings import settings

celery_app = Celery(
    __name__, broker=settings.celery_broker, backend=settings.celery_backend
)
celery_app.conf.enable_utc = True


def send_task(name: str, payload: BaseModel):
    return celery_app.send_task(name=name, args=[payload.model_dump()])


async def wait_async(async_result: AsyncResult) -> Any:
    while not async_result.ready():
        await asyncio.sleep(0.5)
    return async_result.get()
