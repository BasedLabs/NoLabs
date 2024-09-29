__all__ = ["ComponentFlow", "serve_workflow"]

import asyncio
import datetime
import uuid
from abc import ABC
from collections import defaultdict, deque
from logging import Logger
from typing import Generic, List, Optional, Dict, Set, Awaitable

from prefect import State, Task, flow, get_run_logger, task, get_client, serve
from prefect.client.schemas.objects import FlowRun, R, StateType, TaskRun, TERMINAL_STATES
from prefect.context import get_run_context
from prefect.deployments import run_deployment
from prefect.states import Cancelled

from domain.exceptions import NoLabsException, ErrorCodes
from nolabs.application import initialize as initialize_application
from nolabs.domain.models.common import Experiment
from nolabs.infrastructure.settings import settings
from workflow.socketio_events_emitter import (
    emit_component_jobs_event,
    emit_finish_component_event,
    emit_finish_job_event,
    emit_start_component_event,
    emit_start_job_event,
)
from nolabs.domain.models.common import ComponentData, Job
from nolabs.domain.workflow.component import (
    Component,
    ComponentTypeFactory,
    Parameter,
    TInput,
    TOutput,
)


class ComponentFlow(ABC, Generic[TInput, TOutput]):
    component_id: uuid.UUID
    component_name: str

    experiment_id: uuid.UUID

    """Whether input was changed after last execution"""

    job_timeout_seconds: Optional[int] = 1
    component_timeout_seconds: Optional[int] = 10

    logger: Logger

    _job_task_inputs_memo = {}

    def __init__(self, component_id: uuid.UUID, component_name: str, experiment_id: uuid.UUID):
        if self.component_timeout_seconds <= self.job_timeout_seconds:
            raise NoLabsException(ErrorCodes.invalid_workflow_timeouts,
                                  message="Component timeout cannot be less or equal to job timeouts")

        self.component_id = component_id
        self.component_name = component_name
        self.logger = get_run_logger()
        self.experiment_id = experiment_id

        self.execute = flow(
            name=component_name,
            flow_run_name=str(component_id),
            timeout_seconds=self.component_timeout_seconds,
            on_running=[self._on_component_running],
            on_crashed=[self._on_component_finished],
            on_failure=[self._on_component_finished],
            on_cancellation=[self._on_component_finished],
            on_completion=[self._on_component_finished],
        )(self.execute)

        self._job_task = task(
            name=f"{self.component_name}-job",
            task_run_name=str(component_id),
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
        flow_run = ctx.flow_run

        extra = {"component_id": self.component_id, "flow_run_id": flow_run.id}

        data: ComponentData = ComponentData.objects.with_id(self.component_id)

        if (data.flow_run_id
                and flow_run.id != data.flow_run_id
                and await is_workflow_running(data.flow_run_id)):
            return Cancelled(message="Instance of this component is already running")

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


async def is_workflow_running(flow_run_id: uuid.UUID) -> bool:
    async with get_client() as client:
        flow_run = await client.read_flow_run(flow_run_id=flow_run_id)

    state = flow_run.state.type

    return state not in TERMINAL_STATES


@flow(
    name="component_flow",
    flow_run_name="{component_id}"
)
async def component_flow(component_id: uuid.UUID):
    initialize_application()
    data: ComponentData = ComponentData.objects.with_id(component_id)
    component: Component = Component.restore(data)
    component_flow = component.component_flow_type(component_id=data.id,
                                                   component_name=component.name,
                                                   experiment_id=data.experiment.id)
    await component_flow.execute()


@flow(
    name='main_flow',
    flow_run_name='{experiment_id}'
)
async def main_flow(experiment_id: uuid.UUID):
    initialize_application()
    experiment = Experiment.objects.with_id(experiment_id)
    context = get_run_context()
    flow_run = context.flow_run

    if (experiment.flow_run_id
            and flow_run.id != experiment.flow_run_id
            and await is_workflow_running(experiment.flow_run_id)):
        return Cancelled(message="Instance of this workflow is already running")

    components = [Component.restore(data=c) for c in ComponentData.objects(experiment=experiment_id)]

    if not components:
        return

    if len(components) == 1:
        component = components[0]
        component_flow = component.component_flow_type(component_id=component.id,
                                                       component_name=component.name,
                                                       experiment_id=experiment_id)
        await component_flow.execute(return_state=True)
        return

    component_map: Dict[uuid.UUID, Component] = {c.id: c for c in components}

    in_degree: Dict[uuid.UUID, int] = defaultdict(int)
    adj_list: Dict[uuid.UUID, Set[uuid.UUID]] = defaultdict(set)

    for component in components:
        for prev_id in component.previous_component_ids:
            adj_list[prev_id].add(component.id)
            in_degree[component.id] += 1

    queue = deque([c.id for c in components if in_degree[c.id] == 0])

    while queue:
        parallel_group = []

        for _ in range(len(queue)):
            current_id = queue.popleft()
            current_component = component_map[current_id]
            parallel_group.append(current_component)

        group = [
            c.component_flow_type(component_id=c.id,
                                  component_name=c.name,
                                  experiment_id=experiment_id).execute(return_state=True)
            for c in parallel_group
        ]

        await asyncio.gather(*group)

        for component in parallel_group:
            for dependent_id in adj_list[component.id]:
                in_degree[dependent_id] -= 1

                if in_degree[dependent_id] == 0:
                    queue.append(dependent_id)


def serve_workflow():
    main_flow_dep = main_flow.to_deployment(name="workflow", version=str(settings.workflow_version))
    component_flow_dep = component_flow.to_deployment(name="component", version=str(settings.workflow_version))
    serve(main_flow_dep, component_flow_dep)


def start_workflow(experiment_id: uuid.UUID) -> Awaitable[FlowRun]:
    return run_deployment(
        name="main_flow/workflow",
        parameters={
            "experiment_id": experiment_id
        }
    )


def start_component_flow(component_id: uuid.UUID) -> Awaitable[FlowRun]:
    return run_deployment(
        name="component_flow/component",
        parameters={
            "component_id": component_id
        }
    )
