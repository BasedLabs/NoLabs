__all__ = [
    'Workflow',
    'Component',
]

import uuid
from typing import Dict, Any, List
from typing import Optional, Type

from mongoengine import Document, UUIDField, ReferenceField, CASCADE, DictField

from nolabs.application.workflow.component import ComponentTypeFactory, Component, ComponentState
from nolabs.application.workflow.definition import WorkflowDefinition, ComponentTemplate, ComponentDefinition, \
    DefaultDefinition, MappingDefinition
from nolabs.application.workflow.executor import DagExecutor
from nolabs.application.workflow.mappings import map_property
from nolabs.application.workflow.properties import ParameterSchema
from nolabs.domain.models.common import Experiment
from nolabs.exceptions import NoLabsException, ErrorCodes


class WorkflowState(Document):
    id: uuid.UUID = UUIDField(primary_key=True)
    experiment: Experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    definition: Dict[str, Any] = DictField()

    meta = {'collection': 'workflows'}

    @staticmethod
    def create(id: uuid.UUID,
               experiment: Experiment,
               definition: WorkflowDefinition) -> 'WorkflowState':
        return WorkflowState(
            id=id,
            experiment=experiment,
            definition=definition.dict()
        )

    def set_definition(self, definition: WorkflowDefinition):
        self.definition = definition.dict()

    def get_definition(self) -> WorkflowDefinition:
        return WorkflowDefinition(**self.definition)


class Workflow:
    _id: uuid.UUID
    _state: WorkflowState

    def __init__(self, id: uuid.UUID, state: WorkflowState):
        self._id = id
        self._state = state

    @classmethod
    def get(cls, id: uuid.UUID) -> Optional['Workflow']:
        state = WorkflowState.objects.with_id(id)

        if not state:
            return None

        return Workflow(id, state)

    def delete(self):
        self._state.delete()

    @property
    def definition(self) -> WorkflowDefinition:
        definition = WorkflowDefinition(**self._state.definition)
        return definition

    @classmethod
    def all_workflow_ids(cls, experiment_id: uuid.UUID) -> List[uuid.UUID]:
        experiment = Experiment.objects.with_id(experiment_id)
        models: List[WorkflowState] = WorkflowState.objects(experiment=experiment).only('id')
        return [m.id for m in models]

    @classmethod
    def create(cls, experiment_id: uuid.UUID) -> 'Workflow':
        experiment: Experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        component_templates: List[ComponentTemplate] = []
        component_definitions: List[ComponentDefinition] = []

        asd = uuid.uuid4()

        for name, bag in ComponentTypeFactory.enumerate():
            component = bag(id=uuid.uuid4(), experiment_id=experiment_id)

            output_schema = component.output_schema
            output_parameters = {name: map_property(prop, output_schema) for name, prop in
                                 output_schema.properties.items()}

            input_schema = component.input_schema
            input_parameters = {name: map_property(prop, input_schema) for name, prop in
                                input_schema.properties.items()}

            component_templates.append(ComponentTemplate(
                name=component.name,
                input=input_parameters,
                output=output_parameters,
                description=component.description
            ))

            if name == 'test1':
                component_definitions.append(
                    ComponentDefinition(
                        name=name,
                        component_id=asd,
                        defaults=[
                            DefaultDefinition(
                                target_path=['x'],
                                value=5
                            ),
                            DefaultDefinition(
                                target_path=['y'],
                                value=15
                            )
                        ]
                    )
                )
            else:
                component_definitions.append(
                    ComponentDefinition(
                        name=name,
                        component_id=uuid.uuid4(),
                        mappings=[
                            MappingDefinition(
                                source_component_id=asd,
                                source_path=['x'],
                                target_path=['y']
                            ),
                            MappingDefinition(
                                source_component_id=asd,
                                source_path=['y'],
                                target_path=['x']
                            )
                        ]
                    )
                )

        id = uuid.uuid4()

        definition = WorkflowDefinition(
            workflow_id=id,
            error=None,
            component_templates=component_templates,
            components=component_definitions
        )

        state = WorkflowState.create(id=id, experiment=experiment, definition=definition)
        state.save(cascade=True)

        return Workflow(id=id, state=state)

    def _get_component(self, component_id: uuid.UUID) -> Optional[Component]:
        state: ComponentState = ComponentState.objects.with_id(component_id)

        if not state:
            return

        return ComponentTypeFactory.get_type(state.name)(
            id=state.id,
            experiment_id=self._state.experiment.id,
            job_ids=state.job_ids,
            input_schema=ParameterSchema(**state.input_schema),
            output_schema=ParameterSchema(**state.output_schema),
            input_value_dict=state.input_value_dict,
            output_value_dict=state.output_value_dict,
            previous_component_ids=state.previous_component_ids
        )

    def update(self, definition: WorkflowDefinition):
        state = self._state

        components: List[Component] = []

        for wf_component in definition.components:
            component = ComponentTypeFactory.get_type(wf_component.name)(id=wf_component.component_id,
                                                                         experiment_id=state.experiment.id)

            components.append(component)

        # Validate graph
        if self._is_cyclic(graph=components):
            definition.error = 'Workflow must be acyclic'
            definition.valid = False
            self._state.set_definition(definition=definition)
            self._state.save(cascade=True)
            return definition

        schema_valid = True

        for workflow_component in definition.components:
            component = [c for c in components if c.id == workflow_component.component_id][0]

            # Validate mappings
            for mapping in workflow_component.mappings:
                source_components = [c for c in definition.components if
                                     c.component_id == mapping.source_component_id]
                if not source_components:
                    raise NoLabsException(ErrorCodes.component_not_found)

                source_component = [c for c in components if c.id == mapping.source_component_id][0]

                component.add_previous(component_id=source_component.id)

                error = component.try_map_property(
                    component=source_component,
                    path_from=mapping.source_path,
                    target_path=mapping.target_path)

                if error:
                    mapping.error = error.msg
                    schema_valid = False

            # Validate defaults

            for default in workflow_component.defaults:
                error = component.try_set_default(target_path=default.target_path, value=default.value)
                if error:
                    default.error = error.msg
                    schema_valid = False

        definition.valid = schema_valid
        self._state.set_definition(definition=definition)
        self._state.save()

        return definition

    async def start(self):
        experiment = self._state.experiment
        definition = self._state.get_definition()

        if not definition.valid:
            raise NoLabsException(ErrorCodes.invalid_workflow_definition)

        components: Dict[uuid.UUID, Component] = {}
        component_states: Dict[uuid.UUID, ComponentState] = {}

        for workflow_component in definition.components:
            c: Component = self._get_component(workflow_component.component_id)
            if c:
                components[workflow_component.component_id] = c
            else:
                component: Component = ComponentTypeFactory.get_type(workflow_component.name)(
                    id=workflow_component.component_id,
                    experiment_id=self._state.experiment.id)
                state = ComponentState.create(
                    id=component.id,
                    workflow=self._state,
                    component=component,
                    job_ids=[]
                )
                components[workflow_component.component_id] = component
                component_states[workflow_component.component_id] = state

        for workflow_component in definition.components:
            component: Component = components[workflow_component.component_id]

            if not component:
                raise NoLabsException(ErrorCodes.component_not_found)

            # Check that component mappings exist
            for mapping in workflow_component.mappings:
                source_component: Component = components[mapping.source_component_id]

                if source_component.id not in component.previous_component_ids:
                    component.add_previous(component_id=source_component.id)

                component.try_map_property(component=source_component,
                                           path_from=mapping.source_path,
                                           target_path=mapping.target_path)

            for default in workflow_component.defaults:
                component.try_set_default(default.target_path, value=default.value)

        for state in component_states.values():
            state.save()

        executor = DagExecutor()
        await executor.execute(workflow_id=self._state.id, experiment_id=experiment.id,
                               components=list(components.values()))

    def _is_cyclic(self, graph: List[Component]) -> bool:
        visited: Set[WorkflowGraphNode] = set()  # type: ignore
        recursion_stack = set()

        def dfs(vertex: Component):
            if vertex in recursion_stack:
                return True
            if vertex in visited:
                return False

            visited.add(vertex)
            recursion_stack.add(vertex)

            for neighbor_id in vertex.previous_component_ids:
                if dfs([n for n in graph if n.id == neighbor_id][0]):
                    return True

            recursion_stack.remove(vertex)
            return False

        for vertex in graph:
            if vertex not in visited:
                if dfs(vertex):
                    return True

        return False

    @classmethod
    def set_component_types(cls, components: List[Type[Component]]):
        types = {c.name: c for c in components}

        ComponentTypeFactory.set_types(types)
