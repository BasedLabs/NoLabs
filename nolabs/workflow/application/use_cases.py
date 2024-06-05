import uuid
from typing import Optional, List, Dict, Type, Union
from uuid import UUID

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Experiment
from nolabs.workflow.application.api_models import GetComponentStateResponse, GetComponentStateRequest, \
    AllWorkflowSchemasResponse
from nolabs.workflow.component import Component
from nolabs.workflow.models import WorkflowSchemaDbModel, ComponentDbModel
from nolabs.workflow.properties import Property, ParameterSchema, Items
from nolabs.workflow.workflow import Workflow
from nolabs.workflow.workflow_schema import WorkflowSchemaModel, ComponentModel, PropertyModel, ItemsModel


def map_items(i: Items, schema: ParameterSchema) -> ItemsModel:
    ref = None

    if i.ref:
        ref = i.ref

    if i.items and i.items.ref:
        ref = i.items.ref

    definition = None

    if ref:
        definition = schema.defs[schema.get_ref_type_name(ref)]

    items = None

    if i.items:
        items = map_items(i.items, schema) if not isinstance(i.items, list) else [map_items(item, schema) for item in
                                                                                  i.items]

    properties = {name: map_property(prop, schema) for name, prop in (definition.properties if
                                                                      definition.properties else {}).items()} if definition else \
        {name: map_property(prop, schema) for name, prop in (i.properties if
                                                             i.properties else {}).items()}

    return ItemsModel(
        type=i.type,
        properties=properties if not items else None,
        required=i.required,
        description=i.description,
        enum=i.enum,
        const=i.const,
        format=i.format,
        default=i.default,
        example=i.example,
        items=items
    )


def map_property(p: Union[Property, Items], schema: ParameterSchema) -> PropertyModel:
    ref = None

    if p.ref:
        ref = p.ref

    if p.items and p.items.ref:
        ref = p.items.ref

    definition = None

    if ref:
        definition = schema.defs[schema.get_ref_type_name(ref)]

    items = None

    if p.items:
        items = map_items(p.items, schema) if not isinstance(p.items, list) else [map_items(item, schema) for item in p.items]

    properties = {name: map_property(prop, schema) for name, prop in (definition.properties if
                                                                        definition.properties else {}).items()} if definition else \
        {name: map_property(prop, schema) for name, prop in (p.properties if
                                                             p.properties else {}).items()}

    return PropertyModel(
        type=p.type,
        properties=properties if not items else None,
        required=p.required,
        description=p.description,
        enum=p.enum,
        const=p.const,
        format=p.format,
        default=p.default,
        example=p.example,
        title=p.title,
        anyOf=[
            map_property(prop, schema=schema)
            for prop in p.anyOf
        ],
        ref=ref,
        items=items
    )

class DeleteWorkflowSchemaFeature:
    async def handle(self, workflow_id: UUID):
        db_model = WorkflowSchemaDbModel.objects.with_id(workflow_id)

        if db_model:
            db_model.delete()

class AllWorkflowSchemasFeature:
    async def handle(self, experiment_id: UUID) -> AllWorkflowSchemasResponse:
        experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        db_models: List[WorkflowSchemaDbModel] = WorkflowSchemaDbModel.objects(experiment=experiment)
        return AllWorkflowSchemasResponse(
            ids=[m.id for m in db_models]
        )

class CreateWorkflowSchemaFeature:
    available_components: Dict[str, Type[Component]]

    def __init__(self, available_components: Dict[str, Type[Component]]):
        self.available_components = available_components

    async def handle(self, experiment_id: UUID) -> WorkflowSchemaModel:
        experiment: Experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        components_models: List[ComponentModel] = []

        for component_name, component_type in self.available_components.items():
            component = component_type(id=uuid.uuid4(), experiment=experiment)

            output_schema = component.output_schema
            output_parameters = {name: map_property(prop, output_schema) for name, prop in
                                 output_schema.properties.items()}

            components_models.append(ComponentModel(
                name=component_name,
                input={},
                output=output_parameters
            ))

        id = uuid.uuid4()

        workflow_schema = WorkflowSchemaModel(
            workflow_id=id,
            error=None,
            components=components_models,
            workflow_components=[]
        )

        db_model = WorkflowSchemaDbModel.create(
            id=id,
            experiment=experiment,
            value=workflow_schema
        )
        db_model.save(cascade=True)

        return workflow_schema

class GetWorkflowSchemaFeature:
    async def handle(self, workflow_id: UUID) -> Optional[WorkflowSchemaModel]:
        db_model: WorkflowSchemaDbModel = WorkflowSchemaDbModel.objects.with_id(workflow_id)

        if not db_model:
            return None

        return db_model.get_workflow_value()

class SetWorkflowSchemaFeature:
    available_components: Dict[str, Type[Component]]

    def __init__(self, available_components: Dict[str, Type[Component]]):
        self.available_components = available_components

    async def handle(self, workflow_schema: WorkflowSchemaModel) -> WorkflowSchemaModel:
        db_model: WorkflowSchemaDbModel = WorkflowSchemaDbModel.objects.with_id(workflow_schema.workflow_id)

        # Validate components names
        for component in workflow_schema.workflow_components:
            if component.name not in [c_name for c_name in self.available_components.keys()]:
                component.error = f'Component with name "{component.name}" not found'
                workflow_schema.valid = False
                db_model.set_workflow_value(workflow_schema)
                db_model.save()
                return workflow_schema

        components: List[Component] = [
            self.available_components[wf.name](id=wf.component_id, experiment=db_model.experiment) for wf in
            workflow_schema.workflow_components
        ]

        # Validate graph
        if self.is_cyclic(graph=components):
            workflow_schema.error = 'Workflow must be acyclic'
            workflow_schema.valid = False
            db_model.set_workflow_value(workflow_schema)
            db_model.save()
            return workflow_schema

        schema_valid = True

        for workflow_component in workflow_schema.workflow_components:
            component = [c for c in components if c.id == workflow_component.component_id][0]

            # Validate mappings
            for mapping in workflow_component.mappings:
                source_components = [c for c in workflow_schema.workflow_components if
                                     c.component_id == mapping.source_component_id]
                if not source_components:
                    raise NoLabsException(ErrorCodes.component_not_found)

                source_component = [c for c in components if c.id == mapping.source_component_id][0]

                component.add_previous(component=source_component)

                error = component.try_map_property(component=source_component,
                                                   path_from=mapping.source_path,
                                                   path_to=mapping.target_path)

                if error:
                    mapping.error = error.msg
                    schema_valid = False

            # Validate defaults

            for default in workflow_component.defaults:
                error = component.try_set_default(path_to=default.path_to, value=default.value)
                if error:
                    default.error = error.msg
                    schema_valid = False

        workflow_schema.valid = schema_valid
        db_model.set_workflow_value(value=workflow_schema)
        db_model.save(cascade=True)

        return workflow_schema

    def is_cyclic(self, graph: List[Component]) -> bool:
        visited: Set[WorkflowGraphNode] = set()  # type: ignore
        recursion_stack = set()

        def dfs(vertex: Component):
            if vertex in recursion_stack:
                return True
            if vertex in visited:
                return False

            visited.add(vertex)
            recursion_stack.add(vertex)

            for neighbor in vertex.previous:
                if dfs(neighbor):  # type: ignore # TODO Change later
                    return True

            recursion_stack.remove(vertex)
            return False

        for vertex in graph:
            if vertex not in visited:
                if dfs(vertex):
                    return True

        return False

class StartWorkflowFeature:
    available_components: Dict[str, Type[Component]]

    def __init__(self, available_components: Dict[str, Type[Component]]):
        self.available_components = available_components

    async def handle(self, experiment_id: UUID):
        experiment = Experiment.objects.with_id(experiment_id)

        db_model: WorkflowSchemaDbModel = WorkflowSchemaDbModel.objects(experiment=experiment).first()

        workflow_schema = db_model.get_workflow_value()

        if not workflow_schema.valid:
            raise NoLabsException(ErrorCodes.invalid_workflow_schema,
                                  f'Run {SetWorkflowSchemaFeature.__name__} to get schema errors')

        components: List[Component] = []

        component: Component
        for workflow_component in workflow_schema.workflow_components:
            component_db_model: ComponentDbModel = ComponentDbModel.objects.with_id(
                workflow_component.component_id)
            if component_db_model:
                component: Component = self.available_components[workflow_component.name](
                    id=workflow_component.component_id,
                    jobs=component_db_model.jobs,
                    input_parameter_dict=component_db_model.input_parameter_dict,
                    output_parameter_dict=component_db_model.output_parameter_dict,
                    experiment=experiment)
            else:
                id = uuid.uuid4()
                component: Component = self.available_components[workflow_component.name](id=id,
                                                                                          experiment=experiment)
                component_db_model = ComponentDbModel(
                    id=id,
                    workflow=db_model,
                    last_exception='',
                    input_parameter_dict={},
                    output_parameter_dict={},
                    jobs=[]
                )
                component_db_model.save()

            components.append(component)

        for workflow_component in workflow_schema.workflow_components:
            component: Component = [c for c in components if c.id == workflow_component.component_id][0]

            if not component:
                raise NoLabsException(ErrorCodes.component_not_found)

            # Check that component mappings exist
            for mapping in workflow_component.mappings:
                source_component: Component = [c for c in components if c.id == mapping.source_component_id][0]

                component.add_previous(component=source_component)
                component.try_map_property(component=source_component,
                                           path_from=mapping.source_path,
                                           path_to=mapping.target_path)

            for default in workflow_component.defaults:
                component.try_set_default(default.path_to, value=default.value)

        workflow = Workflow()
        await workflow.execute(workflow_schema=workflow_schema, components=components)

class GetComponentParametersFeature:
    async def handle(self, request: GetComponentStateRequest) -> GetComponentStateResponse:
        component: ComponentDbModel = ComponentDbModel.objects.with_id(request.component_id)

        if not component:
            raise NoLabsException(ErrorCodes.component_not_found)

        return GetComponentStateResponse(
            input_dict=component.input_parameter_dict,
            output_dict=component.output_parameter_dict,
            job_ids=[j.iid.value for j in component.jobs]
        )
