import json
import uuid
from abc import abstractmethod, ABC
from typing import Dict, Optional, Any

from celery.result import AsyncResult

from nolabs.domain.exceptions import NoLabsException, ErrorCodes
from infrastructure.celery_app_factory import get_celery_app
from infrastructure.redis_client_factory import Redis
from nolabs.infrastructure.redis_client_factory import use_redis_pipe
from nolabs.workflow.core.states import ControlStates, celery_to_internal_mapping
from nolabs.workflow.core.states import state_transitions

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
        state = await Redis.client.get(cid)
        if not state:
            self._state_cache = ControlStates.UNKNOWN
            return self._state_cache
        self._state_cache = ControlStates(state)
        return self._state_cache

    async def get_input(self) -> Optional[Dict[str, Any]]:
        cid = f"{self._id}:input"
        execution_input = await Redis.client.get(cid)
        if not execution_input:
            return None
        return json.loads(execution_input)

    async def set_input(self, execution_input: Dict[str, Any]):
        cid = f"{self._id}:input"
        await Redis.client.set(cid, json.dumps(execution_input, default=str))

    async def set_output(self, output: Dict[str, Any]):
        cid = f"{self._id}:output"
        await Redis.client.set(cid, output or "")

    async def get_output(self) -> Dict[str, Any]:
        cid = f"{self._id}:output"
        return await Redis.client.get(cid)

    async def set_message(self, message: str):
        cid = f"{self._id}:message"
        await Redis.client.set(cid, message or "")

    async def get_message(self):
        cid = f"{self._id}:message"
        return await Redis.client.get(cid)

    async def set_state(self, state: ControlStates):
        current_state = self._state_cache
        if current_state and state not in state_transitions[current_state]:
            raise NoLabsException(ErrorCodes.invalid_states_transition, f"Cannot move from {current_state} to {state}")

        cid = f"{self._id}:state"
        await Redis.client.set(cid, state)
        self._state_cache = state

    @abstractmethod
    async def sync_started(self):
        ...

    async def can_schedule(self, **kwargs) -> bool:
        """Schedule and start as soon as possible"""
        state = await self.get_state()
        return state == ControlStates.UNKNOWN

    async def schedule(self, **kwargs):
        await self.set_state(state=ControlStates.SCHEDULED)

    async def can_start(self, **kwargs):
        return await self.get_state() == ControlStates.SCHEDULED

    @abstractmethod
    async def start(self, **kwargs):
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
        celery_task_id = await Redis.client.get(cid)
        if not celery_task_id:
            raise NoLabsException(ErrorCodes.celery_task_startd_but_task_id_not_found)

        async_result = AsyncResult(id=celery_task_id, app=celery)
        celery_state = async_result.state
        new_state = celery_to_internal_mapping[celery_state]

        if new_state != state:
            result = async_result.result
            async with use_redis_pipe():
                await self.set_state(state=new_state)
                if new_state in [ControlStates.FAILURE, ControlStates.CANCELLED]:
                    await self.set_message(message=str(result))
                    await self.set_output(output={})
                else:
                    await self.set_output(output=result)

        return new_state

    async def _prepare_for_start(self) -> str:
        """
        :returns: celery_task_id
        """
        cid = self.celery_task_id_cid
        task_id = str(uuid.uuid4())
        await Redis.client.set(cid, task_id)
        return task_id