# Convention: service.method


import asyncio
import uuid
from typing import Any

from celery import Celery
from celery.result import AsyncResult, allow_join_result
from pydantic import BaseModel

import microservices.esmfold_light.service.api_models as esmfold_light
from nolabs.infrastructure.settings import settings


class Cel:
    _app: Celery

    def __init__(self):
        self._app = Celery(
            __name__, broker=settings.celery_broker, backend=settings.celery_backend
        )
        self._app.conf.update(
            enable_utc=True,
            task_default_queue='workflow',
            task_track_started=True,
            worker_send_task_events=True,
            task_acks_late=True,
            task_reject_on_worker_lost=True,
            worker_state_db="/opt/celery-state.db",
            accept_content=["application/json", "application/x-python-serialize"],
            task_always_eager=False,
            task_eager_propagates=False
        )
        self._app.autodiscover_tasks(force=True)

    @property
    def app(self) -> Celery:
        return self._app

    def _send_task(self, name: str, payload: BaseModel) -> AsyncResult:
        return self._app.send_task(name=name, args=[payload.model_dump()])

    def send_task(
        self, id: uuid.UUID, name: str, queue: str, payload: BaseModel
    ) -> AsyncResult:
        return self._app.send_task(
            id=str(id), name=name, queue=queue, args=[payload.model_dump()]
        )

    async def _wait_async(self, async_result: AsyncResult) -> Any:
        while not async_result.ready():
            await asyncio.sleep(0.5)
        with allow_join_result():
            return async_result.get()

    def cancel_task(self, task_id: str):
        self._app.control.revoke(task_id)

    async def esmfold_light_inference(
        self, task_id: str, payload: esmfold_light.InferenceInput
    ) -> esmfold_light.InferenceOutput:
        async_result = self._app.send_task(
            id=task_id,
            name="esmfold-light-service.inference",
            queue="esmfold-light-service",
            args=[payload.model_dump()],
        )
        result = await self._wait_async(async_result)
        return esmfold_light.InferenceOutput(**result)

cel = Cel()
