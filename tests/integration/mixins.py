import asyncio
import time
import uuid
from typing import Type, Dict, Tuple, List, Optional

from celery.result import AsyncResult

from nolabs.domain.models.common import Experiment, ExperimentName, ExperimentId, ComponentData
from nolabs.infrastructure.celery_app_factory import get_celery_app
from nolabs.workflow.core.component import Component, ComponentTypeFactory
from nolabs.workflow.core.graph import Graph
from nolabs.workflow.core.states import TERMINAL_STATES


class SeedExperimentMixin:
    async def seed_experiment(self, id: uuid.UUID):
        e = Experiment.create(id=ExperimentId(id), name=ExperimentName("Test experiment"))
        e.save()


class SeedComponentsMixin:
    def seed_component(self, experiment_id: uuid.UUID, component_type: Type[Component], component_id: Optional[uuid.UUID] = None) -> Component:
        ComponentTypeFactory.add_type(component_type)
        component = component_type(id=component_id or uuid.uuid4())
        data = ComponentData.create(id=component.id, experiment=experiment_id)
        component.dump(data=data)
        data.save()
        return component

    def seed_mappings(self, component: Component, previous_components: List[Tuple[Component, List[str], List[str]]]):
        for prev_c, map_from, map_to in previous_components:
            component.add_previous(prev_c.id)
            component.try_map_property(prev_c, map_from, map_to)
        data = ComponentData.objects.with_id(component.id)
        component.dump(data=data)
        data.save()

class GraphTestMixin:
    async def await_for_celery_task(self, task_id: str):
        celery = get_celery_app()
        async_result = AsyncResult(id=task_id, app=celery)
        while True:
            ready = async_result.ready()
            if ready:
                return
            await asyncio.sleep(0.1)
