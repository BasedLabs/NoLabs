# Convention: service.method


import asyncio
import uuid
from typing import Any

import microservices.diffdock.service.api_models as diffdock
import microservices.esmfold_light.service.api_models as esmfold_light
import microservices.reinvent.service.api_models as reinvent
from celery import Celery
from celery.result import AsyncResult, allow_join_result
from pydantic import BaseModel

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
            worker_state_db="celery-state.db",
            accept_content=["application/json", "application/x-python-serialize"]
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
        self, task_id: str, request: reinvent.RunSamplingRequest
    ):
        task = self._app.send_task(
            id=task_id,
            name="reinvent.run_sampling",
            queue="reinvent",
            args=[request.model_dump()],
        )
        await self._wait_async(task)

    async def reinvent_run_learning(
        self, task_id: str, request: reinvent.RunReinforcementLearningRequest
    ):
        task = self._app.send_task(
            id=task_id,
            name="reinvent.run_reinforcement_learning",
            queue="reinvent",
            args=[request.model_dump()],
        )
        await self._wait_async(task)

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
        self, task_id: str, payload: diffdock.RunDiffDockPredictionRequest
    ) -> diffdock.RunDiffDockPredictionResponse:
        async_result = self._app.send_task(
            id=task_id,
            name="diffdock.inference",
            queue="diffdock",
            args=[payload.model_dump()],
        )
        result = await self._wait_async(async_result)
        return diffdock.RunDiffDockPredictionResponse(**result)

    async def reinvent_prepare_target(
        self, task_id: str, payload: reinvent.PreparePdbqtRequest
    ) -> reinvent.PreparePdbqtResponse:
        task = self._app.send_task(
            id=task_id,
            name="reinvent.prepare_target",
            queue="reinvent",
            args=[payload.model_dump()],
        )
        result = await self._wait_async(task)
        return reinvent.PreparePdbqtResponse(**result)


cel = Cel()
