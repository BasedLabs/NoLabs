import json
import uuid
from abc import abstractmethod, ABC
from typing import List, Set, Dict, Optional, Any

from celery.result import AsyncResult

from nolabs.infrastructure.cel import cel
from nolabs.infrastructure.log import get_scheduler_logger
from nolabs.infrastructure.red import redis_client, use_redis_pipe
from nolabs.workflow.component import Component
from nolabs.workflow.logic.states import ControlStates, celery_to_internal_mapping
from workflow.logic.component_execution_nodes import ComponentExecutionNode

logger = get_scheduler_logger()


class ExecutionNode(ABC):
    """
    Represents single distributed task on a worker
    """

    def __init__(self, id: str):
        self._id = id

    async def get_state(self) -> Optional[ControlStates]:
        cid = f"{self._id}:state"
        state = await redis_client.get(cid)
        if not state:
            return ControlStates.UNKNOWN
        return ControlStates(state)

    async def get_input(self) -> Optional[Dict[str, Any]]:
        cid = f"{self._id}:execution_input"
        execution_input = await redis_client.get(cid)
        if not execution_input:
            return None
        return json.loads(execution_input)

    async def set_input(self, execution_input: Dict[str, Any]):
        cid = f"{self._id}:input"
        await redis_client.set(cid, json.dumps(execution_input))

    async def set_output(self, output: Dict[str, Any]):
        cid = f"{self._id}:output"
        await redis_client.set(cid, output)

    async def get_output(self) -> Dict[str, Any]:
        cid = f"{self._id}:output"
        return await redis_client.get(cid)

    async def set_error(self, error: str):
        cid = f"{self._id}:error"
        await redis_client.set(cid, error)

    async def get_error(self):
        cid = f"{self._id}:error"
        await redis_client.get(cid)

    async def set_state(self, state: ControlStates):
        cid = f"{self._id}:state"
        await redis_client.set(cid, state)

    @abstractmethod
    async def synchronize(self):
        ...

    @abstractmethod
    async def can_execute(self) -> bool:
        ...

    @abstractmethod
    async def execute(self, **kwargs):
        ...

    async def schedule(self, **kwargs):
        await self.set_state(state=ControlStates.SCHEDULED)


class CeleryExecutionNode(ExecutionNode, ABC):
    async def synchronize(self) -> ControlStates:
        state = await self.get_state()
        if state != ControlStates.STARTED:
            return state
        cid = f"{self._id}:celery_task_id"
        async with use_redis_pipe():
            celery_task_id = await redis_client.get(cid)
            if not celery_task_id:
                logger.critical("Node in trackable state, however celery_task_id is empty")
                await self.set_state(state=ControlStates.FAILURE)
                await self.set_error(error="Internal error")
                return ControlStates.FAILURE
            async_result = AsyncResult(id=celery_task_id, app=cel.app)
            new_state = celery_to_internal_mapping[async_result.state]

            if new_state != state:
                await self.set_state(state=new_state)
                if new_state in [ControlStates.FAILURE, ControlStates.CANCELLED]:
                    await self.set_error(error=async_result.result)
                await self.set_output(output=async_result.result)

    async def _prepare_for_execute(self) -> str:
        """
        :returns: celery_task_id
        """
        async with use_redis_pipe():
            cid = f"{self._id}:state"
            redis_client.set(cid, ControlStates.STARTED)
            cid = f"{self._id}:celery_task_id"
            task_id = str(uuid.uuid4())
            redis_client.set(cid, task_id)
            return task_id


class ExecutionGraph(ExecutionNode):
    def __init__(self, experiment_id: uuid.UUID, id: str):
        super().__init__(f"graph:{experiment_id}")

    async def synchronize(self):
        pass

    async def execute(self, **kwargs):
        graph = await self.get_input()
        for component_id in graph.keys():
            node = ComponentExecutionNode(uuid.UUID(component_id))
            await node.execute()

    async def can_execute(self) -> bool:
        return await self.get_state() == ControlStates.SCHEDULED

    async def schedule(self, components: List[Component]):
        graph = {}

        for component in components:
            graph[str(component.id)] = set()
            for component_2 in components:
                if component_2.id == component.id:
                    continue
                if component_2.id in component.previous_component_ids:
                    graph[str(component.id)].add(component_2.id)

        async with use_redis_pipe():
            await self.set_input(execution_input=graph)
            await self.set_state(state=ControlStates.SCHEDULED)


class ExecutionGraph:
    def __init__(self, experiment_id: uuid.UUID, graph: Dict[uuid.UUID, Set[uuid.UUID]]):
        self._experiment_id = experiment_id
        self._graph = graph

    @classmethod
    async def graph_exists(cls, experiment_id: uuid.UUID) -> bool:
        return await redis_client.exists(f"graph:{experiment_id}")

    @classmethod
    async def get(cls, experiment_id: uuid.UUID) -> Optional['ExecutionGraph']:
        instance = await redis_client.get(f"graph:{experiment_id}")
        if instance:
            return ExecutionGraph(experiment_id, json.loads(instance))
        return None

    @classmethod
    async def build(cls, experiment_id: uuid.UUID, components: List[Component]) -> 'ExecutionGraph':
        graph = {}

        for component in components:
            graph[component.id] = set()
            for component_2 in components:
                if component_2.id == component.id:
                    continue
                if component_2.id in component.previous_component_ids:
                    graph[component.id].add(component_2.id)

        await redis_client.set(f"graph:{experiment_id}", json.dumps(graph))
        return ExecutionGraph(experiment_id, graph)

    async def start(self):
        pass

    async def stop(self):
        pass



    async def monitor(self):
        pass

    def __getitem__(self, component_id: uuid.UUID) -> Set[uuid.UUID]:
        if component_id not in self._graph:
            return set()

        return self._graph[component_id]
