import asyncio
import time
import uuid
from typing import Type

from nolabs.domain.models.common import Experiment, ExperimentName, ExperimentId, ComponentData
from nolabs.workflow.core.component import Component, ComponentTypeFactory
from nolabs.workflow.core.graph import GraphExecutionNode
from workflow.core.states import TERMINAL_STATES


class SeedExperimentMixin:
    async def seed_experiment(self, id: uuid.UUID):
        e = Experiment.create(id=ExperimentId(id), name=ExperimentName("Test experiment"))
        e.save()


class SeedComponentsMixin:
    def seed_component(self, experiment_id: uuid.UUID, component_type: Type[Component]) -> Component:
        ComponentTypeFactory.add_type(component_type)
        component = component_type(id=uuid.uuid4())
        data = ComponentData.create(id=component.id, experiment=experiment_id)
        component.dump(data=data)
        data.save()
        return component

    def seed_direct_components(self, experiment_id: uuid.UUID, from_c: Component, to_c: Component):
        pass

class GraphTestMixin:
    async def sync_until_terminal(self, graph: GraphExecutionNode, timeout=20):
        t = time.time()
        while await graph.get_state() not in TERMINAL_STATES and time.time() - t < timeout:
            await graph.sync_started()
            await asyncio.sleep(0.1)

