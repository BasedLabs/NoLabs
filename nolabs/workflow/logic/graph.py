import json
import uuid
from abc import abstractmethod, ABC
from typing import List, Set, Dict, Optional, Any

from celery.result import AsyncResult

from infrastructure.cel import cel
from infrastructure.log import get_scheduler_logger
from infrastructure.red import redis_client, use_redis_pipe
from workflow.component import Component
from workflow.logic.control import _component_main_task, _complete_component_task
from workflow.logic.job_execution_nodes import JobMainTaskExecutionNode
from workflow.logic.states import ControlStates, celery_to_internal_mapping


NON_PROPAGATE_STATES = [ControlStates.FAILURE, ControlStates.CANCELLED]
PROPAGATE_STATES = [ControlStates.SUCCESS]
PROGRESS_STATES = [ControlStates.STARTED]

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
        pipe = redis_client.pipeline()
        cid = f"{self._id}:state"
        pipe.set(cid, ControlStates.STARTED)
        cid = f"{self._id}:celery_task_id"
        task_id = str(uuid.uuid4())
        pipe.set(cid, task_id)
        await pipe.execute()
        return task_id


class ComponentExecutionNode(ExecutionNode):
    def __init__(self, experiment_id: uuid.UUID, component_id: uuid.UUID):
        super().__init__(id=f"graph:{component_id}")
        self.experiment_id = experiment_id
        self.component_id = component_id

    async def synchronize(self):
        main_task = ComponentMainTaskExecutionNode(experiment_id=self.experiment_id,
                                                   component_id=self.component_id)
        await main_task.synchronize()
        complete_component_task = ComponentCompleteTaskExecutionNode(experiment_id=self.experiment_id,
                                                                     component_id=self.component_id)
        await complete_component_task.synchronize()

    async def get_state(self) -> ControlStates:
        steps: List[ExecutionNode] = [
            ComponentMainTaskExecutionNode(experiment_id=self.experiment_id,
                                           component_id=self.component_id),
            ComponentCompleteTaskExecutionNode(experiment_id=self.experiment_id,
                                               component_id=self.component_id)
        ]

        states = []
        for step in steps:
            states.append(await step.get_state())

        if ControlStates.FAILURE in states:
            return ControlStates.FAILURE

    async def can_execute(self) -> bool:
        main_task = ComponentMainTaskExecutionNode(experiment_id=self.experiment_id, component_id=self.component_id)

        job_ids = await main_task.get_jobs()
        for job_id in job_ids:
            job_task = JobMainTaskExecutionNode

        complete_task = ComponentCompleteTaskExecutionNode(experiment_id=self.experiment_id,
                                                           component_id=self.component_id)

        return await main_task.can_execute() or await complete_task.can_execute()

    async def move(self):
        main_task = ComponentMainTaskExecutionNode(experiment_id=self.experiment_id, component_id=self.component_id)
        if await main_task.should_move():
            await main_task.move()
            return

        complete_task = ComponentCompleteTaskExecutionNode(experiment_id=self.experiment_id,
                                                           component_id=self.component_id)
        if await complete_task.should_move():
            await complete_task.move()

    async def allow_next(self):
        complete_task = ComponentCompleteTaskExecutionNode(experiment_id=self.experiment_id,
                                                           component_id=self.component_id)
        return await complete_task.allow_next()

    async def next(self) -> List['ExecutionNode']:
        graph = await ExecutionGraph.get(experiment_id=self.experiment_id)
        next_items = []
        for item in graph[self.component_id]:
            next_items.append(ComponentExecutionNode(experiment_id=self.experiment_id, component_id=item))
        return next_items


class ComponentMainTaskExecutionNode(CeleryExecutionNode):
    def __init__(self, experiment_id: uuid.UUID, component_id: uuid.UUID):
        super().__init__(id=f"graph:{component_id}:{_component_main_task.name}")
        self.component_id = component_id
        self.experiment_id = experiment_id

    async def should_move(self) -> bool:
        state = await self.get_state()

        if state in TERMINAL_STATES:
            return False

        if state == ControlStates.STARTED:
            return False

        previous_component_ids = await ExecutionGraph.get(experiment_id=self.experiment_id)
        for prev_comp_id in previous_component_ids:
            node = ComponentExecutionNode(experiment_id=self.experiment_id, component_id=prev_comp_id)
            if not await node.allow_propagate():
                return False
        return True

    async def move(self):
        celery_task_id = await self._prepare_for_start()
        _component_main_task.apply_async(
            kwargs={"experiment_id": self.experiment_id, "component_id": self.component_id},
            task_id=celery_task_id,
            retry=False)

    async def allow_propagate(self) -> bool:
        ...

    async def get_jobs(self) -> List[uuid.UUID]:
        return ....


class ComponentCompleteTaskExecutionNode(CeleryExecutionNode):
    def __init__(self, experiment_id: uuid.UUID, component_id: uuid.UUID):
        super().__init__(id=f"graph:{component_id}:{_complete_component_task.name}")
        self.component_id = component_id
        self.experiment_id = experiment_id

    async def should_execute(self) -> bool:
        if await self.is_final():
            return False

        if await self.started():
            return False

        main_task = ComponentMainTaskExecutionNode(
            experiment_id=self.experiment_id,
            component_id=self.component_id
        )

        if not await main_task.is_final() or not await main_task.allow_propagate():
            return False

        for job_id in await main_task.get_jobs():
            job_task = JobMainTaskExecutionNode(
                experiment_id=self.experiment_id,
                component_id=self.component_id,
                job_id=job_id
            )

            if not await job_task.is_final() or not await job_task.allow_propagate():
                return False

        return True

    async def execute(self):
        main_task = ComponentMainTaskExecutionNode(
            experiment_id=self.experiment_id,
            component_id=self.component_id
        )

        job_ids = await main_task.get_jobs()
        celery_task_id = await self._prepare_for_start()
        _complete_component_task.apply_async(
            kwargs={"experiment_id": self.experiment_id, "component_id": self.component_id, "job_ids": job_ids},
            task_id=celery_task_id,
            retry=False)


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

    def __getitem__(self, component_id: uuid.UUID) -> Set[uuid.UUID]:
        if component_id not in self._graph:
            return set()

        return self._graph[component_id]
