import asyncio
from functools import lru_cache
from typing import Any

from celery import Celery
from celery.result import AsyncResult

from nolabs.infrastructure.settings import settings


@lru_cache
def get_celery_app() -> Celery:
    app = Celery(
        __name__, broker=settings.celery_broker, backend=settings.celery_backend
    )
    app.conf.update(
        timezone = 'Europe/London',
        enable_utc=True,
        worker_pool=settings.celery_worker_pool,
        task_default_queue="workflow",
        task_track_started=True,
        worker_send_task_events=True,
        task_acks_late=True,
        task_reject_on_worker_lost=True,
        # worker_state_db=settings.celery_worker_state_db,
        accept_content=["application/json", "application/x-python-serialize"],
        task_always_eager=False,
        task_eager_propagates=False,
        redbeat_redis_url=settings.celery_backend
    )
    app.autodiscover_tasks(force=True)

    return app

async def wait_for_task(task_id: str) -> Any:
    celery = get_celery_app()
    async_result = AsyncResult(id=task_id, app=celery)
    while True:
        ready = async_result.ready()
        if ready:
            break
        await asyncio.sleep(0.1)
    return async_result.get()
