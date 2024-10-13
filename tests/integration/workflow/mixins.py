import time
import uuid
from typing import Type

from domain.models.common import Experiment, ExperimentName, ExperimentId, ComponentData
from workflow.core.component import Component
from workflow.core.graph import GraphExecutionNode


class SeedExperimentMixin:
    async def seed_experiment(self, id: uuid.UUID):
        e = Experiment.create(id=ExperimentId(id), name=ExperimentName("Test experiment"))
        e.save()


class SeedComponentsMixin:
    def seed_component(self, experiment_id: uuid.UUID, component_type: Type[Component]) :
        component = component_type(id=uuid.uuid4())
        data = ComponentData.create(id=component.id, experiment=experiment_id)
        component.dump(data=data)
        data.save()
        return component

    def seed_direct_components(self, experiment_id: uuid.UUID, from_c: Component, to_c: Component):
        pass

class GraphTestMixin:
    async def sync_until_terminal(self, graph: GraphExecutionNode, timeout=10):
        t = time.time()
        while time.time() - t < timeout:
            await graph.sync_started()

