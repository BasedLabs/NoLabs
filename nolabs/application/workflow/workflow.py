__all__ = [
    'create_schema',
    'update_schema',
    'delete',
    'get_schema'
]

import uuid
from typing import Dict, List, Any
from typing import Optional, Type

from pydantic import BaseModel

from nolabs.application.workflow.component import ComponentTypeFactory, Component, Parameter
from nolabs.application.workflow.data import WorkflowData, ComponentData
from nolabs.application.workflow.mappings import map_property
from nolabs.application.workflow.prefect.dag import PrefectDagExecutor
from nolabs.application.workflow.schema import WorkflowSchema, ComponentSchemaTemplate, ComponentSchema
from nolabs.exceptions import NoLabsException, ErrorCodes


# TODO make it functional


def get_component(component_id: uuid.UUID) -> Optional[Component]:
    state: ComponentData = ComponentData.objects.with_id(component_id)

    if not state:
        return

    return ComponentTypeFactory.get_type(state.name)(
        id=state.id,
        input_schema=Parameter(**state.input_schema),
        output_schema=Parameter(**state.output_schema),
        input_value_dict=state.input_value_dict,
        output_value_dict=state.output_value_dict,
        previous_component_ids=state.previous_component_ids
    )


def is_job_running(job_id: uuid.UUID) -> bool:
    jobs = ComponentData.objects(jobs_runs__running=True).only('jobs_runs')


class InputPropertyErrorView(BaseModel):
    loc: List[str]
    msg: str


class JobInputErrorView(BaseModel):
    job_id: uuid.UUID
    msg: str


class JobExceptionView(BaseModel):
    job_id: uuid.UUID
    msg: str


class ComponentRunView(BaseModel):
    input_dict: Dict[str, Any]
    output_dict: Dict[str, Any]
    job_ids: List[uuid.UUID]
    input_property_errors: List[InputPropertyErrorView]
    exception: str
    jobs_input_errors: List[JobInputErrorView]
    jobs_exception: List[JobExceptionView]


def get_latest_component_run(component_id: uuid.UUID) -> Optional[ComponentRunView]:
    state = ComponentData.objects.with_id(component_id)

    if not state:
        return None


def create_schema(id: uuid.UUID) -> 'WorkflowSchema':
    component_templates: List[ComponentSchemaTemplate] = []
    component_schemas: List[ComponentSchema] = []

    for name, bag in ComponentTypeFactory.enumerate():
        component = bag(id=uuid.uuid4())

        output_schema = component.output_schema
        output_parameters = {name: map_property(prop, output_schema) for name, prop in
                             output_schema.properties.items()}

        input_schema = component.input_schema
        input_parameters = {name: map_property(prop, input_schema) for name, prop in
                            input_schema.properties.items()}

        component_templates.append(ComponentSchemaTemplate(
            name=component.name,
            input=input_parameters,
            output=output_parameters,
            description=component.description
        ))

    schema = WorkflowSchema(
        workflow_id=id,
        error=None,
        component_templates=component_templates,
        components=component_schemas
    )

    state = WorkflowData.create(id=id, schema=schema)
    state.save(cascade=True)

    return schema


def get_schema(id: uuid.UUID) -> WorkflowSchema:
    data: WorkflowData = WorkflowData.objects.with_id(id)
    return data.get_schema()


def delete_schema(id: uuid.UUID):
    WorkflowData.objects.with_id(id).delete()


def set_component_types(components: List[Type[Component]]):
    types = {c.name: c for c in components}

    ComponentTypeFactory.set_types(types)


def update_schema(schema: WorkflowSchema) -> WorkflowSchema:
    data = WorkflowData.objects.get(id=schema.workflow_id)

    components: List[Component] = []

    for wf_component in schema.components:
        component = ComponentTypeFactory.get_type(wf_component.name)(id=wf_component.component_id)

        components.append(component)

    # Validate graph
    if _is_cyclic(graph=components):
        schema.error = 'Workflow must be acyclic'
        schema.valid = False
        data.set_schema(schema=schema)
        data.save(cascade=True)
        return schema

    schema_valid = True

    for workflow_component in schema.components:
        component = [c for c in components if c.id == workflow_component.component_id][0]

        # Validate mappings
        for mapping in workflow_component.mappings:
            source_components = [c for c in schema.components if
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

    schema.valid = schema_valid
    data.set_schema(schema=schema)
    data.save()

    return schema


def _is_cyclic(graph: List[Component]) -> bool:
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


async def start(id: uuid.UUID, extra_dag_parameters: Optional[Dict[str, Any]] = None):
    data = WorkflowData.objects.with_id(id)

    schema = data.get_schema()

    if not schema.valid:
        raise NoLabsException(ErrorCodes.invalid_workflow_schema)

    components: List[Component] = [Component.restore(ComponentData.objects.with_id(c.component_id)) for c in
                                   schema.components]
    executor = PrefectDagExecutor()
    await executor.execute(workflow_id=id,
                           components=components,
                           extra=extra_dag_parameters)
