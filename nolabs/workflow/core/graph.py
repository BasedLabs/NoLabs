import asyncio
import time
import uuid
from typing import List, Dict, Optional

from pydantic import BaseModel
from redis.client import Pipeline

from nolabs.domain.exceptions import NoLabsException, ErrorCodes
from nolabs.infrastructure.celery_app_factory import get_celery_app
from nolabs.infrastructure.log import logger
from nolabs.infrastructure.redis_client_factory import redlock, get_redis_pipe
from nolabs.workflow.core import Tasks
from nolabs.workflow.core.component import Component
from nolabs.workflow.core.component_execution_nodes import ComponentExecutionNode
from nolabs.workflow.core.job_execution_nodes import JobExecutionNode
from nolabs.workflow.core.node import ExecutionNode
from nolabs.workflow.core.states import ControlStates, PROGRESS_STATES
from nolabs.workflow.core.states import ERROR_STATES
from nolabs.workflow.core.states import TERMINAL_STATES


class GraphData(BaseModel):
    graph: Dict[uuid.UUID, List[uuid.UUID]]
    schedule: List[uuid.UUID]


class Graph(ExecutionNode):
    def __init__(self, experiment_id: uuid.UUID):
        super().__init__(f"execution_node:{experiment_id}")
        self.experiment_id = experiment_id

    async def sync_started(self):
        state = await self.get_state()

        if state != ControlStates.STARTED:
            return

        model = GraphData(**await self.get_input())

        all_completed = True

        for component_id in model.schedule:
            node = ComponentExecutionNode(experiment_id=self.experiment_id, component_id=component_id)
            state = await node.get_state()

            previous_component_ids = model.graph[component_id]

            should_cancel = False

            # Cancel if node is not running and not terminal and any previous component is failed
            if state != ControlStates.STARTED and state not in TERMINAL_STATES:
                for previous_component_id in previous_component_ids:
                    prev_node = ComponentExecutionNode(experiment_id=self.experiment_id,
                                                       component_id=previous_component_id)
                    prev_node_state = await prev_node.get_state()
                    if prev_node_state in ERROR_STATES:
                        should_cancel = True

                if should_cancel:
                    await node.set_cancelled(reason="Previous component failed")
                    continue

            if await node.can_start(previous_component_ids=previous_component_ids):
                await node.start()

            await node.sync_started()
            node_state = await node.get_state()
            # By checking for Unknown state we skipped not scheduled components
            if node_state not in TERMINAL_STATES and node_state != ControlStates.UNKNOWN:
                all_completed = False

        if all_completed:
            await self.set_state(ControlStates.SUCCESS)

    async def start(self, **kwargs):
        lock = redlock(key=self._id, blocking=True)
        await lock.acquire()

        try:
            if not await self.can_start():
                raise NoLabsException(ErrorCodes.start_workflow_failed)

            await self.set_state(ControlStates.STARTED)
        finally:
            if await lock.locked():
                await lock.release()

    async def set_components_graph(self, components: List[Component]):
        # TODO ensure that is not running
        lock = redlock(key=self._id, blocking=True)
        await lock.acquire()

        try:
            state = await self.get_state()

            if state in PROGRESS_STATES:
                raise NoLabsException(ErrorCodes.workflow_running)

            graph = {}

            for component in components:
                graph[component.id] = []
                for component_2 in components:
                    if component_2.id == component.id:
                        continue
                    if component_2.id in component.previous_component_ids:
                        graph[component.id].append(component_2.id)

            model = GraphData(
                graph=graph,
                schedule=[]
            )
            await self.set_input(execution_input=model)
        finally:
            if await lock.locked():
                await lock.release()

    async def can_schedule(self) -> bool:
        state = await self.get_state()
        return state != ControlStates.CANCELLING

    async def schedule(self, schedule: List[Component]):
        lock = redlock(key=self._id, blocking=True)
        await lock.acquire()

        try:
            if not await self.can_schedule():
                raise NoLabsException(ErrorCodes.cannot_schedule_node)

            inp = await self.get_input()
            model = GraphData(**inp)
            model.schedule = list(set(model.schedule + [c.id for c in schedule]))

            if not set(model.schedule).issubset(set(model.graph.keys())):
                raise NoLabsException(ErrorCodes.scheduled_components_are_not_in_graph)

            for component in schedule:
                node = ComponentExecutionNode(experiment_id=self.experiment_id, component_id=component.id)
                if await node.can_schedule():
                    await node.schedule()

            await self.set_input(execution_input=model)
            state = await self.get_state()

            if state == ControlStates.UNKNOWN:
                await super().schedule()
                return

            if state in [ControlStates.STARTED, ControlStates.SCHEDULED]:
                return

            if state in TERMINAL_STATES:
                pipe = get_redis_pipe()
                await super().reset(pipe=pipe)
                await super().schedule()
        finally:
            if await lock.locked():
                await lock.release()

    def get_component_node(self, component_id: uuid.UUID) -> ComponentExecutionNode:
        return ComponentExecutionNode(experiment_id=self.experiment_id, component_id=component_id)

    def get_job_node(self, component_id: uuid.UUID, job_id: uuid.UUID) -> JobExecutionNode:
        return JobExecutionNode(experiment_id=self.experiment_id, component_id=component_id, job_id=job_id)

    async def sync_cancelling(self):
        if await self.get_state() != ControlStates.CANCELLING:
            return

        inp = await self.get_input()
        model = GraphData(**inp)

        all_cancelled = True

        for component_id in model.graph:
            node = ComponentExecutionNode(experiment_id=self.experiment_id, component_id=component_id)
            if await node.can_cancel():
                await node.cancel()
            node_state = await node.sync_cancelling()
            if node_state not in TERMINAL_STATES:
                all_cancelled = False

        if all_cancelled:
            await self.set_state(ControlStates.CANCELLED)

    async def reset(self, pipe: Pipeline):
        inp = await self.get_input()
        model = GraphData(**inp)

        pipe = get_redis_pipe()
        for component_id in model.graph:
            node = ComponentExecutionNode(experiment_id=self.experiment_id, component_id=component_id)
            await node.reset(pipe=pipe)
        await super().reset(pipe=pipe)
        await pipe.execute()

    async def sync(self):
        logger.info('Graph sync check', extra={
            'experiment_id': self.experiment_id
        })
        lock = redlock(key=self._id, blocking=False, auto_release_time=10.0)

        if not await lock.acquire():
            logger.info('Graph is already syncing, returning', extra={
                'experiment_id': self.experiment_id
            })
            return

        try:
            state = await self.get_state()

            if state in TERMINAL_STATES or state == ControlStates.UNKNOWN:
                return

            logger.info('Graph is syncing', extra={
                'experiment_id': self.experiment_id
            })

            if await self.can_start():
                await self.start()

            if await self.get_state() in ControlStates.STARTED:
                await self.sync_started()

            if await self.get_state() in ControlStates.CANCELLING:
                await self.sync_cancelling()
        except Exception as e:
            pipe = get_redis_pipe()
            await self.set_state(state=ControlStates.FAILURE, pipe=pipe)
            await self.set_message(message=str(e), pipe=pipe)
            await pipe.execute()
            logger.exception()
        finally:
            logger.info('Graph syncing finish', extra={
                'experiment_id': self.experiment_id
            })
            if await lock.locked():
                await lock.release()
