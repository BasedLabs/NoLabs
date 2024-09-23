__all__ = ["ComponentFlow"]

import datetime
import uuid
from abc import ABC
from logging import Logger
from typing import Generic, List, Optional

from prefect import State, Task, flow, get_run_logger, task
from prefect.cache_policies import INPUTS
from prefect.client.schemas.objects import FlowRun, R, StateType, TaskRun
from prefect.context import get_run_context
from prefect.states import Completed

from nolabs.application.workflow.api.socketio_events_emitter import (
    emit_component_jobs_event, emit_finish_component_event,
    emit_finish_job_event, emit_start_component_event, emit_start_job_event)
from nolabs.application.workflow.component import (Component,
                                                   ComponentTypeFactory,
                                                   Parameter, TInput, TOutput)
from nolabs.domain.models.common import ComponentData, Job


def _name_builder(name: str, id: uuid.UUID):
    return f"{name},{str(id)}"


def _run_name_builder(name: str, at: datetime.datetime):
    return f"{name},At:{at.isoformat()}"


class ComponentFlow(ABC, Generic[TInput, TOutput]):
    component_id: uuid.UUID
    component_name: str

    experiment_id: uuid.UUID

    """Whether input was changed after last execution"""

    job_timeout_seconds: Optional[int] = 1
    component_timeout_seconds: Optional[int] = 10

    logger: Logger

    _job_task_inputs_memo = {}

    def __init__(self, component: Component, experiment_id: uuid.UUID):
        self.component_id = component.id
        self.component_name = component.name

        self.logger = get_run_logger()

        execute_flow_name = _name_builder(
            id=self.component_id, name=self.component_name
        )

        self.experiment_id = experiment_id

        self.execute = flow(
            name=execute_flow_name,
            flow_run_name=_run_name_builder(
                name=execute_flow_name, at=datetime.datetime.utcnow()
            ),
            timeout_seconds=self.component_timeout_seconds,
            on_running=[self._on_component_running],
            on_crashed=[self._on_component_finished],
            on_failure=[self._on_component_finished],
            on_cancellation=[self._on_component_finished],
            on_completion=[self._on_component_finished],
        )(self.execute)

        execute_task_name = _name_builder(
            id=self.component_id, name=f"{self.component_name}-job"
        )

        execute_task_run_name = execute_task_name + ",{job_id},{at}"

        self._job_task = task(
            name=execute_task_name,
            task_run_name=execute_task_run_name,
            timeout_seconds=self.job_timeout_seconds,
            on_failure=[self._on_task_finish],
            on_completion=[self._on_task_finish],
        )(self._job_task)

    async def get_jobs(self, inp: TInput) -> List[uuid.UUID]:
        return []

    async def gather_jobs(
        self, inp: TInput, job_ids: List[uuid.UUID]
    ) -> Optional[TOutput]:
        pass

    def _on_task_finish(self, task: Task, task_run: TaskRun, state: State):
        emit_finish_job_event(
            experiment_id=self.experiment_id,
            component_id=self.component_id,
            job_id=self._job_task_inputs_memo[task_run.id],
        )

    async def _job_task(self, job_id: uuid.UUID, at: datetime.datetime) -> State[R]:
        ctx = get_run_context()
        task_run_id = ctx.task_run.id

        self._job_task_inputs_memo[task_run_id] = job_id

        emit_start_job_event(
            experiment_id=self.experiment_id,
            component_id=self.component_id,
            job_id=job_id,
        )

        job: Job = Job.objects.with_id(job_id)
        job.set_run_info(task_run_id=task_run_id, executed_at=at)
        await job.save()
        return await self.job_task(job_id=job_id)

    async def job_task(self, job_id: uuid.UUID) -> State[R]:
        pass

    async def execute(self):
        ctx = get_run_context()
        flow_run_id = ctx.flow_run.id

        extra = {"component_id": self.component_id, "flow_run_id": flow_run_id}

        data: ComponentData = ComponentData.objects.with_id(self.component_id)
        component = Component.restore(data=data)

        prev_components: List[Component] = []

        for previous_component_id in data.previous_component_ids:
            previous_component_state = ComponentData.objects.with_id(
                previous_component_id
            )
            previous_component = self._get_component(
                from_state=previous_component_state
            )

            errors = previous_component.output_errors()
            if errors:
                self.logger.info(
                    "Previous component output errors",
                    extra={
                        **extra,
                        **{
                            "previous_component_id": previous_component_id,
                            "errors": [(e.msg, e.loc) for e in errors],
                        },
                    },
                )
                return

            prev_components.append(previous_component)

        input_changed = component.set_input_from_previous(prev_components)

        if input_changed:
            self.logger.info("Input changed, checking input", extra=extra)

            input_errors = component.input_errors()

            if input_errors:
                self.logger.info(
                    "Input errors",
                    extra={
                        **extra,
                        **{"input_errors": [(e.msg, e.loc) for e in input_errors]},
                    },
                )

                return

        input_value = component.input_value
        job_ids = await self.get_jobs(inp=input_value)

        extra = {**extra, **{"job_ids": job_ids}}

        self.logger.info("Retrieved job ids", extra=extra)

        job_errors = False

        if job_ids:
            emit_component_jobs_event(
                experiment_id=self.experiment_id,
                component_id=self.component_id,
                job_ids=job_ids,
            )

            states = self._job_task.map(
                job_ids, at=datetime.datetime.utcnow(), return_state=True
            )
            job_errors = any([s for s in states if s.type != StateType.COMPLETED])

            if job_errors:
                self.logger.info("Jobs errors", extra=extra)

        if not job_errors:
            output = await self.gather_jobs(inp=input_value, job_ids=job_ids)
            component.output_value = output or {}
        component.dump(data=data)
        data.save()

    def _get_component(self, from_state: ComponentData) -> Component:
        component = ComponentTypeFactory.get_type(from_state.name)(
            id=from_state.id,
            input_schema=Parameter(**from_state.input_schema),
            output_schema=Parameter(**from_state.output_schema),
            input_value_dict=from_state.input_value_dict,
            output_value_dict=from_state.output_value_dict,
            previous_component_ids=from_state.previous_component_ids,
        )

        return component

    def _on_component_running(self, _, flow_run: FlowRun, state: State):
        data: ComponentData = ComponentData.objects.with_id(self.component_id)
        data.flow_run_id = flow_run.id
        data.last_executed_at = datetime.datetime.utcnow()
        data.prefect_state = state.type
        data.save()

        emit_start_component_event(self.experiment_id, self.component_id)

    def _on_component_finished(self, _, flow_run: FlowRun, state: State):
        emit_finish_component_event(
            experiment_id=self.experiment_id, component_id=self.component_id
        )
