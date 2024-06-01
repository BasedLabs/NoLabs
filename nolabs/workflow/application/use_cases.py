import uuid
from typing import Optional, List, Dict, Type
from uuid import UUID

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Experiment
from nolabs.workflow.application.api_models import GetComponentStateResponse, GetComponentStateRequest, \
    AllWorkflowSchemasResponse
from nolabs.workflow.component import PythonComponent
from nolabs.workflow.models import WorkflowSchemaDbModel, PythonComponentDbModel
from nolabs.workflow.properties import Property
from nolabs.workflow.workflow import Workflow
from nolabs.workflow.workflow_schema import WorkflowSchemaModel, ComponentModel, PropertyModel


def map_property(property: Property) -> PropertyModel:
    return PropertyModel(
        type=property.type,
        properties={name: map_property(prop) for name, prop in (property.properties if
                                                                property.properties else {}).items()},
        required=property.required,
        description=property.description,
        enum=property.enum,
        const=property.const,
        format=property.format,
        default=property.default,
        example=property.example,
        title=property.title,
        anyOf=[
            map_property(prop)
            for prop in property.anyOf
        ],
        ref=property.ref
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
    available_components: Dict[str, Type[PythonComponent]]

    def __init__(self, available_components: Dict[str, Type[PythonComponent]]):
        self.available_components = available_components

    async def handle(self, experiment_id: UUID) -> WorkflowSchemaModel:
        experiment: Experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        components_models: List[ComponentModel] = []

        for component_name, component_type in self.available_components.items():
            component = component_type(id=uuid.uuid4(), experiment=experiment)
            components_models.append(ComponentModel(
                name=component_name,
                input={name: map_property(prop) for name, prop in component.input_properties.items()},
                output={name: map_property(prop) for name, prop in component.output_properties.items()}
            ))

        workflow_schema = WorkflowSchemaModel(
            workflow_id=uuid.uuid4(),
            error=None,
            components=components_models,
            workflow_components=[]
        )

        db_model = WorkflowSchemaDbModel.create(
            id=uuid.uuid4(),
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
    available_components: Dict[str, Type[PythonComponent]]

    def __init__(self, available_components: Dict[str, Type[PythonComponent]]):
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

        components: List[PythonComponent] = [
            PythonComponent(id=wf.component_id, experiment=db_model.experiment) for wf in
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

    def is_cyclic(self, graph: List[PythonComponent]) -> bool:
        visited: Set[WorkflowGraphNode] = set()  # type: ignore
        recursion_stack = set()

        def dfs(vertex: PythonComponent):
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
    available_components: Dict[str, Type[PythonComponent]]

    def __init__(self, available_components: Dict[str, Type[PythonComponent]]):
        self.available_components = available_components

    async def handle(self, experiment_id: UUID):
        experiment = Experiment.objects.with_id(experiment_id)

        db_model: WorkflowSchemaDbModel = WorkflowSchemaDbModel.objects(experiment=experiment).first()

        workflow_schema = db_model.get_workflow_value()

        if not workflow_schema.valid:
            raise NoLabsException(ErrorCodes.invalid_workflow_schema,
                                  f'Run {SetWorkflowSchemaFeature.__name__} to get schema errors')

        components: List[PythonComponent] = []

        component: PythonComponent
        for workflow_component in workflow_schema.workflow_components:
            component_db_model: PythonComponentDbModel = PythonComponentDbModel.objects.with_id(
                workflow_component.component_id)
            if component_db_model:
                component: PythonComponent = PythonComponent(id=workflow_component.component_id,
                                                             jobs=component_db_model.jobs,
                                                             input_parameter_dict=component_db_model.input_parameter_dict,
                                                             output_parameter_dict=component_db_model.output_parameter_dict,
                                                             experiment=experiment)
            else:
                id = uuid.uuid4()
                component: PythonComponent = self.available_components[workflow_component.name](id=id,
                                                                                                experiment=experiment)
                component_db_model = PythonComponentDbModel(
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
            component: PythonComponent = [c for c in components if c.id == workflow_component.component_id][0]

            if not component:
                raise NoLabsException(ErrorCodes.component_not_found)

            # Check that component mappings exist
            for mapping in workflow_component.mappings:
                source_component: PythonComponent = [c for c in components if c.id == mapping.source_component_id][0]

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
        component: PythonComponentDbModel = PythonComponentDbModel.objects.with_id(request.component_id)

        if not component:
            raise NoLabsException(ErrorCodes.component_not_found)

        return GetComponentStateResponse(
            input_dict=component.input_parameter_dict,
            output_dict=component.output_parameter_dict,
            job_ids=[j.iid.value for j in component.jobs]
        )