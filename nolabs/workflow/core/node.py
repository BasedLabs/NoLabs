import json
import uuid
from abc import abstractmethod, ABC
from typing import Dict, Optional, Any

from celery.result import AsyncResult
from redis.asyncio import Redis

from domain.exceptions import NoLabsException, ErrorCodes
from infrastructure.celery_app_factory import get_celery_app
from infrastructure.redis_client_factory import get_redis_client
from nolabs.infrastructure.log import get_worker_logger
from nolabs.infrastructure.redis_client_factory import use_redis_pipe
from nolabs.workflow.core.states import ControlStates, celery_to_internal_mapping
from nolabs.workflow.core.states import state_transitions

logger = get_worker_logger()
redis_client = get_redis_client()
celery = get_celery_app()

class ExecutionNode(ABC):
    """
    Represents single distributed task on a worker
    """

    def __init__(self, id: str):
        self._id = id
        self._state_cache = None

    async def get_state(self) -> Optional[ControlStates]:
        if self._state_cache:
            return self._state_cache

        cid = f"{self._id}:state"
        state = await redis_client.get(cid)
        if not state:
            self._state_cache = ControlStates.UNKNOWN
            return self._state_cache
        self._state_cache = ControlStates(state)
        return self._state_cache

    async def get_input(self) -> Optional[Dict[str, Any]]:
        cid = f"{self._id}:execution_input"
        execution_input = await redis_client.get(cid)
        if not execution_input:
            return None
        return json.loads(execution_input)

    async def set_input(self, execution_input: Dict[str, Any]):
        cid = f"{self._id}:input"
        await redis_client.set(cid, json.dumps(execution_input, default=str))

    async def set_output(self, output: Dict[str, Any]):
        cid = f"{self._id}:output"
        await redis_client.set(cid, output or "")

    async def get_output(self) -> Dict[str, Any]:
        cid = f"{self._id}:output"
        return await redis_client.get(cid)

    async def set_message(self, message: str):
        cid = f"{self._id}:message"
        await redis_client.set(cid, message or "")

    async def get_message(self):
        cid = f"{self._id}:error"
        return await redis_client.get(cid)

    async def set_state(self, state: ControlStates):
        current_state = await self.get_state()
        if state not in state_transitions[current_state]:
            raise NoLabsException(ErrorCodes.invalid_states_transition, f"Cannot move from {current_state} to {state}")

        cid = f"{self._id}:state"
        await redis_client.set(cid, state)
        self._state_cache = state

    @abstractmethod
    async def sync_started(self):
        ...

    async def can_schedule(self, **kwargs) -> bool:
        """Schedule and execute as soon as possible"""
        return await self.get_state() == ControlStates.UNKNOWN

    async def schedule(self, **kwargs):
        await self.set_state(state=ControlStates.SCHEDULED)

    async def can_execute(self, **kwargs):
        return await self.get_state() == ControlStates.SCHEDULED

    @abstractmethod
    async def execute(self, **kwargs):
        ...

    async def started(self) -> bool:
        return await self.get_state() == ControlStates.STARTED


class CeleryExecutionNode(ExecutionNode, ABC):
    def __init__(self, id: str):
        super().__init__(id)

    @property
    def celery_task_id_cid(self) -> str:
        return f"{self._id}:celery_task_id"

    async def sync_started(self) -> ControlStates:
        state = await self.get_state()
        if state != ControlStates.STARTED:
            return state
        cid = self.celery_task_id_cid
        async with use_redis_pipe():
            celery_task_id = await redis_client.get(cid)
            if not celery_task_id:
                raise NoLabsException(ErrorCodes.celery_task_executed_but_task_id_not_found)
            async_result = AsyncResult(id=celery_task_id, app=celery)
            celery_state = async_result.state
            new_state = celery_to_internal_mapping[celery_state]

            if new_state != state:
                result = async_result.get()
                await self.set_state(state=new_state)
                if new_state in [ControlStates.FAILURE, ControlStates.CANCELLED]:
                    await self.set_message(message=result)
                await self.set_output(output=result)

    async def _prepare_for_execute(self) -> str:
        """
        :returns: celery_task_id
        """
        async with use_redis_pipe():
            cid = self.celery_task_id_cid
            task_id = str(uuid.uuid4())
            await redis_client.set(cid, task_id)
            return task_id