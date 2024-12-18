import json
import uuid
from typing import Dict, List, Optional

from pydantic import BaseModel

from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.infrastructure.log import logger
from nolabs.infrastructure.redis_client_factory import rd, redlock
from nolabs.workflow.core.component import Component
from nolabs.workflow.core.component_execution_nodes import ComponentExecutionNode
from nolabs.workflow.core.job_execution_nodes import JobExecutionNode
from nolabs.workflow.core.states import (
    ERROR_STATES,
    PROGRESS_STATES,
    TERMINAL_STATES,
    ControlStates,
)


class GraphMetadata(BaseModel):
    graph: Dict[uuid.UUID, List[uuid.UUID]]
    schedule: List[uuid.UUID]


class Graph:
    def __init__(self, experiment_id: uuid.UUID):
        self._id = f"execution_node:{experiment_id}"
        self.experiment_id = experiment_id

    def get_component_node(self, component_id: uuid.UUID) -> ComponentExecutionNode:
        return ComponentExecutionNode(
            experiment_id=self.experiment_id, component_id=component_id
        )

    def get_job_node(
        self, component_id: uuid.UUID, job_id: uuid.UUID
    ) -> JobExecutionNode:
        return JobExecutionNode(
            experiment_id=self.experiment_id, component_id=component_id, job_id=job_id
        )

    async def started(self) -> bool:
        metadata = await self._get_metadata()
        if not metadata:
            return False
        for component_id in metadata.graph:
            node = self.get_component_node(component_id)
            if await node.get_state() in PROGRESS_STATES:
                return True

        return False

    async def set_components_graph(self, components: List[Component]):
        lock = redlock(key=self._id, blocking=True)
        lock.acquire()

        try:
            metadata = await self._get_metadata()

            if metadata:
                for component_id in metadata.graph:
                    node = self.get_component_node(component_id)
                    if await node.get_state() in PROGRESS_STATES:
                        raise NoLabsException(ErrorCodes.workflow_running)

            graph = {}

            for component in components:
                graph[component.id] = []
                for component_2 in components:
                    if component_2.id == component.id:
                        continue
                    if component_2.id in component.previous_component_ids:
                        graph[component.id].append(component_2.id)

            metadata = GraphMetadata(graph=graph, schedule=[])
            await self._set_metadata(data=metadata)
        finally:
            if lock.owned():
                lock.release()

    async def schedule(self, component_ids: List[uuid.UUID]):
        lock = redlock(key=self._id, blocking=True)

        try:
            metadata = await self._get_metadata()
            metadata.schedule = list(set(metadata.schedule + component_ids))

            if not set(metadata.schedule).issubset(set(metadata.graph.keys())):
                raise NoLabsException(ErrorCodes.scheduled_components_are_not_in_graph)

            for component_id in component_ids:
                node = ComponentExecutionNode(
                    experiment_id=self.experiment_id, component_id=component_id
                )
                if await node.can_schedule():
                    await node.schedule()

            await self._set_metadata(data=metadata)
        finally:
            if lock.owned():
                lock.release()

    async def _set_metadata(self, data: GraphMetadata):
        cid = f"{self._id}:input"
        v = data.model_dump_json()
        rd.set(cid, v)

    async def metadata_empty(self):
        cid = f"{self._id}:input"
        return rd.exists(cid)

    async def _get_metadata(self) -> Optional[GraphMetadata]:
        cid = f"{self._id}:input"
        metadata = rd.get(cid)
        if not metadata:
            return None
        return GraphMetadata(**json.loads(metadata))

    async def cancel(self):
        metadata = await self._get_metadata()
        if not metadata:
            return
        for component_id in metadata.graph:
            node = self.get_component_node(component_id)
            if await node.can_cancel():
                await node.cancel()

    async def sync(self):
        lock = redlock(key=self._id, blocking=False, auto_release_time=10.0)

        if not lock.acquire():
            logger.debug(
                "Graph is already syncing, returning",
                extra={"experiment_id": self.experiment_id},
            )
            return

        try:
            metadata = await self._get_metadata()
            if not metadata:
                return

            for component_id in metadata.graph:
                node = ComponentExecutionNode(
                    experiment_id=self.experiment_id, component_id=component_id
                )

                state = await node.get_state()

                if state in TERMINAL_STATES:
                    continue

                if state == ControlStates.CANCELLING:
                    await node.sync_cancelling()
                    continue

                previous_component_ids = metadata.graph[component_id]

                should_cancel = False

                # Cancel if node is not running and not terminal and any previous component is failed
                if state != ControlStates.STARTED:
                    for previous_component_id in previous_component_ids:
                        prev_node = self.get_component_node(previous_component_id)
                        prev_node_state = await prev_node.get_state()
                        if prev_node_state in ERROR_STATES:
                            should_cancel = True

                    if should_cancel:
                        await node.set_cancelled(reason="Previous component failed")
                        continue

                if await node.can_start(previous_component_ids=previous_component_ids):
                    await node.start()

                await node.sync_started()
        finally:
            if lock.owned():
                lock.release()
