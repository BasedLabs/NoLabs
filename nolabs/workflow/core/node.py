import json
import uuid
from abc import abstractmethod, ABC
from typing import Dict, Optional, Any

from celery.result import AsyncResult
from pydantic import BaseModel
from redis.client import Pipeline

from nolabs.infrastructure.redis_client_factory import get_redis_pipe
from nolabs.domain.exceptions import NoLabsException, ErrorCodes
from nolabs.infrastructure.celery_app_factory import get_celery_app
from nolabs.infrastructure.redis_client_factory import Redis
from nolabs.workflow.core.states import ControlStates, celery_to_internal_mapping
from nolabs.workflow.core.states import state_transitions
from nolabs.workflow.core.states import ERROR_STATES, TERMINAL_STATES


class ExecutionNode(ABC):
    """
    Represents single distributed task on a worker
    """

    def __init__(self, id: str):
        self._id = id
        self._state_cache = None
        self._state_key = f"{self._id}:state"

    async def get_state(self) -> Optional[ControlStates]:
        #if self._state_cache:
        #    return self._state_cache

        state = await Redis.client.get(self._state_key)
        if not state:
            self._state_cache = ControlStates.UNKNOWN
            return self._state_cache
        self._state_cache = ControlStates(state)
        return self._state_cache

    async def on_started(self):
        ...

    async def on_finished(self):
        ...

    async def get_input(self) -> Optional[Dict[str, Any]]:
        cid = f"{self._id}:input"
        execution_input = await Redis.client.get(cid)
        if not execution_input:
            return None
        return json.loads(execution_input)

    async def set_input(self, execution_input: Dict[str, Any] | str | BaseModel, pipe:Optional[Pipeline] = None):
        cid = f"{self._id}:input"

        data = ''

        if isinstance(execution_input, dict):
            data = json.dumps(execution_input, default=str)
        if isinstance(execution_input, BaseModel):
            data = execution_input.model_dump_json()
        if isinstance(execution_input, str):
            data = execution_input

        await (pipe or Redis.client).set(cid, data)

    async def set_output(self, output: Dict[str, Any], pipe:Optional[Pipeline] = None):
        cid = f"{self._id}:output"
        await (pipe or Redis.client).set(cid, output or "")

    async def get_output(self) -> Dict[str, Any]:
        cid = f"{self._id}:output"
        return await Redis.client.get(cid)

    async def set_message(self, message: str, pipe:Optional[Pipeline] = None):
        cid = f"{self._id}:message"
        await (pipe or Redis.client).set(cid, message or "")

    async def get_message(self):
        cid = f"{self._id}:message"
        return await Redis.client.get(cid)

    async def set_cancelled(self, reason: Optional[str] = None):
        pipe = get_redis_pipe()
        await self.set_state(state=ControlStates.CANCELLED, pipe=pipe)
        await self.set_message(message=reason or '', pipe=pipe)
        await pipe.execute()

    async def set_state(self, state: ControlStates, pipe:Optional[Pipeline] = None):
        current_state = self._state_cache
        if current_state and state not in state_transitions[current_state]:
            raise NoLabsException(ErrorCodes.invalid_states_transition, f"Cannot move from {current_state} to {state}")

        await (pipe or Redis.client).set(self._state_key, state)
        self._state_cache = state

        if state == ControlStates.STARTED:
            await self.on_started()

        if state in TERMINAL_STATES:
            await self.on_finished()

    @abstractmethod
    async def sync_started(self):
        ...

    @abstractmethod
    async def sync_cancelling(self):
        ...

    async def can_schedule(self, **kwargs) -> bool:
        """Schedule and start as soon as possible"""
        state = await self.get_state()
        return state == ControlStates.UNKNOWN

    async def schedule(self, **kwargs):
        await self.set_state(state=ControlStates.SCHEDULED, pipe=(kwargs['pipe'] if 'pipe' in kwargs else None))

    async def can_start(self, **kwargs):
        return await self.get_state() == ControlStates.SCHEDULED

    @abstractmethod
    async def start(self, **kwargs):
        ...

    async def started(self) -> bool:
        return await self.get_state() == ControlStates.STARTED

    async def can_reset(self) -> bool:
        state = await self.get_state()
        return state not in [ControlStates.STARTED, ControlStates.CANCELLED]

    async def reset(self, pipe: Pipeline):
        await self.set_state(ControlStates.UNKNOWN, pipe=pipe)
        await self.set_output(output={}, pipe=pipe)
        await self.set_message(message="", pipe=pipe)

    async def cancel(self):
        if await self.can_cancel():
            await self.set_state(state=ControlStates.CANCELLING)

    async def can_cancel(self):
        state = await self.get_state()
        return state == ControlStates.STARTED


class CeleryExecutionNode(ExecutionNode, ABC):
    def __init__(self, id: str):
        super().__init__(id)
        self.celery = get_celery_app()

    @property
    def celery_task_id_cid(self) -> str:
        return f"{self._id}:celery_task_id"

    async def sync_cancelling(self):
        state = await self.get_state()
        if state != ControlStates.CANCELLING:
            return state
        cid = self.celery_task_id_cid
        celery_task_id = await Redis.client.get(cid)
        if not celery_task_id:
            pipe = get_redis_pipe()
            await self.set_state(state=ControlStates.CANCELLED, pipe=pipe)
            await self.set_output(output={}, pipe=pipe)
            await self.set_message(message="", pipe=pipe)
            return

        celery = get_celery_app()
        celery.control.revoke(celery_task_id, terminate=True)

        async_result = AsyncResult(id=celery_task_id, app=self.celery)
        celery_state = async_result.state
        new_state = celery_to_internal_mapping[celery_state]

        if new_state in TERMINAL_STATES:
            pipe = get_redis_pipe()
            await self.set_state(state=ControlStates.CANCELLED, pipe=pipe)
            await self.set_output(output={}, pipe=pipe)
            await self.set_message(message="", pipe=pipe)
            await pipe.execute()

    async def sync_started(self) -> ControlStates:
        state = await self.get_state()
        if state != ControlStates.STARTED:
            return state
        cid = self.celery_task_id_cid
        celery_task_id = await Redis.client.get(cid)
        if not celery_task_id:
            raise NoLabsException(ErrorCodes.celery_task_executed_but_task_id_not_found)

        async_result = AsyncResult(id=celery_task_id, app=self.celery)
        celery_state = async_result.state
        new_state = celery_to_internal_mapping[celery_state]

        if new_state != state:
            result = async_result.result
            pipe = get_redis_pipe()
            await self.set_state(state=new_state, pipe=pipe)
            if new_state in ERROR_STATES:
                await self.set_message(message=str(result),pipe=pipe)
                await self.set_output(output={},pipe=pipe)
            else:
                await self.set_output(output=result,pipe=pipe)
            await pipe.execute()

        return new_state

    async def _prepare_for_start(self, pipe:Optional[Pipeline] = None) -> str:
        """
        :returns: celery_task_id
        """
        cid = self.celery_task_id_cid
        task_id = str(uuid.uuid4())
        await (pipe or Redis.client).set(cid, task_id)
        return task_id