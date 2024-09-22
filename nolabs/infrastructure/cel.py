# Convention: service.method


import asyncio
import uuid
from typing import Any, Optional, Tuple

import microservices.diffdock.service.api_models as diffdock
import microservices.esmfold_light.service.api_models as esmfold_light
import microservices.reinvent.service.api_models as reinvent
from celery import Celery
from celery.result import AsyncResult
from pydantic import BaseModel

from nolabs.infrastructure.settings import settings


class Cel:
    _app: Celery

    def __init__(self):
        self._app = Celery(
            __name__, broker=settings.celery_broker, backend=settings.celery_backend
        )
        self._app.conf.enable_utc = True
        self._app.autodiscover_tasks(force=True)

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
        return async_result.get()

    async def wait_async(self, async_result: AsyncResult) -> Any:
        while not async_result.ready():
            await asyncio.sleep(0.5)
        return async_result.get()

    async def task_result(self, task_id: str) -> AsyncResult:
        return AsyncResult(id=task_id, app=self._app)

    def cancel_task(self, task_id: str):
        self._app.control.revoke(task_id)

    @property
    def ready_states(self) -> Any:
        return self._app.backend.READY_STATES

    @property
    def failed_states(self) -> Any:
        return self._app.backend.FAILED_STATES

    async def reinvent_run_sampling(
        self, request: reinvent.RunSamplingRequest, wait=True
    ) -> Tuple[Any, str]:
        task = self._send_task("reinvent.run_sampling", payload=request)
        id = task.id
        if not wait:
            return (None, id)
        result = await self._wait_async(task)
        return (result, id)

    async def reinvent_run_learning(
        self, request: reinvent.RunReinforcementLearningRequest, wait=True
    ) -> Tuple[Any, str]:
        task = self._send_task("reinvent.run_reinforcement_learning", payload=request)
        id = task.id
        if not wait:
            return (None, id)
        result = await self._wait_async(task)
        return (result, id)

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

    async def diffdock_inference(
        self, payload: diffdock.RunDiffDockPredictionRequest, wait=True
    ) -> Tuple[Optional[diffdock.RunDiffDockPredictionResponse], str]:
        task = self._send_task("diffdock.inference", payload=payload)
        id = task.id
        if not wait:
            return (None, id)
        result = await self._wait_async(task)
        return (diffdock.RunDiffDockPredictionResponse(**result), id)

    async def reinvent_prepare_target(
        self, payload: reinvent.PreparePdbqtRequest, wait=True
    ) -> Tuple[Optional[reinvent.PreparePdbqtResponse], str]:
        task = self._send_task("reinvent.prepare_target", payload=payload)
        id = task.id
        if not wait:
            return (None, id)
        result = await self._wait_async(task)
        return (reinvent.PreparePdbqtResponse(**result), id)


cel = Cel()
