import uuid
from abc import ABC
from typing import Generic, List, Optional

from celery import Celery
from celery.result import AsyncResult
from celery.states import FAILURE

from nolabs.domain.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.common import ComponentData
from nolabs.domain.workflow.component import TInput, TOutput, Component
from nolabs.infrastructure.log import logger
from nolabs.infrastructure.settings import settings

# tasks

# start -> spin ups overriden start, that will
# gather() -> TOutput

_app = Celery(
    __name__, broker=settings.celery_broker, backend=settings.celery_backend
)
_app.conf.update(
    enable_utc=True,
    task_track_started=True,
    task_acks_late=True,
    task_reject_on_worker_lost=True,
    worker_state_db="celery-state.db",
    accept_content=["application/json", "application/x-python-serialize"]
)
_app.autodiscover_tasks(force=True)


@_app.task(name="celery_flows.setup_component_task_id_task")
def setup_component_task_id_task(component_id: uuid.UUID,
                       task_id: uuid.UUID):
    data: ComponentData = ComponentData.objects.with_id(component_id)

    if data.task_id is not None:
        async_result = AsyncResult(id=data.task_id, app=_app)
        if not async_result.ready():
            raise NoLabsException(ErrorCodes.component_running)

    data.set_task_id(task_id=task_id)
    data.save()

@_app.task(name="celery_flows.on_unexpected_exception")
def on_component_unexpected_exception(request, exc: Exception, traceback: str, component_id: uuid.UUID):
    logger.exception(exc)
    data: ComponentData = ComponentData.objects.with_id(component_id)
    data.set_state(state=FAILURE, state_message=str(exc))
    data.save()

@_app.task(time_limit=600, name="celery_flows.start_component_task")
def start_component_task(component_id: uuid.UUID, experiment_id: uuid.UUID):
    data: ComponentData = ComponentData.objects.with_id(component_id)

    if data.task_id is not None:
        async_result = AsyncResult(id=data.task_id, app=_app)
        if not async_result.ready():
            raise NoLabsException(ErrorCodes.component_running)

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


def start_component_chain()


class _ComponentFlow(Generic[TInput, TOutput]):
    def __init__(self,
                 component_id: uuid.UUID,
                 experiment_id: uuid.UUID):
        self.component_id = component_id
        self.experiment_id = experiment_id

    # will be called from within task
    def started(self, inp: TInput):
        pass

    # call this to complete the comp
    def complete(self,
                 output: Optional[TOutput] = None,
                 message: Optional[str] = None):
        data: ComponentData = ComponentData.objects.with_id(self.component_id)

    # call this to fail the comp
    def fail(self, e: Exception, message: str):
        pass

    # call this to cancel the comp
    def cancel(self, message: str):
        pass


class CeleryComponentLifecycle(_ComponentLifecycle[TInput, TOutput]):
    def setup_jobs(self) -> List[uuid.UUID]:
        pass

    def gather_jobs(self, job_ids: List[uuid.UUID]):
        pass


class _JobLifecycle:
    def __init__(self, job_id: uuid.UUID):
        self.job_id = job_id

    # will be called from within the task
    def start(self):
        pass

    # call this to complete the comp
    def completed(self, message: str):
        pass

    # call this to fail the comp
    def failed(self, e: Exception, message: str):
        pass

    # call this to cancel the comp
    def cancelled(self, message: str):
        pass


class LocalJobLifecycle(_JobLifecycle):
    pass


class CeleryJobLifecycle(_JobLifecycle):
    pass
