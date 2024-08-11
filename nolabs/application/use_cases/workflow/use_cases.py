import logging
import uuid
from typing import Optional, List, Dict, Type
from uuid import UUID

from nolabs.application.use_cases.workflow.api_models import (GetComponentStateResponse, GetComponentStateRequest,
                                                              AllWorkflowDefinitionsResponse, JobErrorResponse,
                                                              InputPropertyErrorResponse, ResetWorkflowRequest,
                                                              StartWorkflowComponentRequest)
from nolabs.application.use_cases.workflow.mappings import map_property
from nolabs.application.workflow.component import Component
from nolabs.application.workflow.definition import (WorkflowDefinition, ComponentTemplate)
from nolabs.application.workflow.models import WorkflowDbModel, ComponentDbModel
from nolabs.domain.models.common import Experiment, Job
from nolabs.exceptions import NoLabsException, ErrorCodes


class DeleteWorkflowDefinitionFeature:
    async def handle(self, workflow_id: UUID):
        db_model = WorkflowDbModel.objects.with_id(workflow_id)

        if db_model:
            db_model.delete()


class AllWorkflowDefinitionsFeature:
    async def handle(self, experiment_id: UUID) -> AllWorkflowDefinitionsResponse:
        experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        ids = WorkflowDbModel.objects(experiment=experiment).only('id')
        return AllWorkflowDefinitionsResponse(
            ids=ids
        )


class CreateWorkflowDefinitionFeature:
    available_components: Dict[str, Type[Component]]

    def __init__(self, available_components: Dict[str, Type[Component]]):
        self.available_components = available_components

    async def handle(self, experiment_id: UUID) -> WorkflowDefinition:
        experiment: Experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        component_templates: List[ComponentTemplate] = []

        for name, bag in self.available_components.items():
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

        id = uuid.uuid4()

        definition = WorkflowDefinition(
            workflow_id=id,
            error=None,
            component_templates=component_templates,
            components=[]
        )

        db_model = WorkflowDbModel.create(
            id=id,
            experiment=experiment,
            components=[],
            definition=definition
        )
        db_model.save(cascade=True)

        return definition


class GetWorkflowDefinitionFeature:
    async def handle(self, workflow_id: UUID) -> Optional[WorkflowDefinition]:
        db_model: WorkflowDbModel = WorkflowDbModel.objects.with_id(workflow_id)

        if not db_model:
            return None

        return db_model.get_workflow_definition()


class UpdateWorkflowSchemaFeature:
    available_components: Dict[str, Type[Component]]

    def __init__(self,
                 available_components: Dict[str, Type[Component]]):
        self.available_components = available_components

    async def handle(self, definition: WorkflowDefinition) -> WorkflowDefinition:
        workflow = self.fetch_workflow(definition.workflow_id)

        components: List[Component] = []

        for wf_component in definition.components:
            component = self.available_components[wf_component.name](id=wf_component.component_id,
                                                                     experiment_id=workflow.experiment.id)

            components.append(component)

        # Validate graph
        if self.is_cyclic(graph=components):
            definition.error = 'Workflow must be acyclic'
            definition.valid = False
            workflow.set_workflow_definition(definition)
            workflow.save()
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
        workflow.set_workflow_definition(value=definition)
        workflow.save(cascade=True)

        return definition

    def fetch_workflow(self, workflow_id) -> WorkflowDbModel:
        return WorkflowDbModel.objects.with_id(workflow_id)

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


class StartWorkflowFeature:
    available_components: Dict[str, Type[Component]]
    logger: logging.Logger

    def __init__(self, available_components: Dict[str, Type[Component]],
                 logger: logging.Logger = logging.getLogger('nolabs')):
        self.available_components = available_components
        self.logger = logger

    async def handle(self, workflow_id: UUID):
        db_model: WorkflowDbModel = self.fetch_workflow(workflow_id)
        experiment = db_model.experiment

        definition = db_model.get_workflow_definition()

        if not definition.valid:
            raise NoLabsException(ErrorCodes.invalid_workflow_definition,
                                  f'Run {UpdateWorkflowSchemaFeature.__name__} to get schema errors')

        components: List[Component] = []

        component: Component
        for workflow_component in definition.components:
            component_db_model: ComponentDbModel = self.fetch_component(workflow_component.component_id)
            if component_db_model:
                component: Component = component_db_model.get_component()
            else:
                id = workflow_component.component_id
                component: Component = self.available_components[workflow_component.name](id=id,
                                                                                                        experiment=experiment)
                component_db_model = ComponentDbModel.create(
                    id=id,
                    workflow=db_model,
                    component=component,
                    jobs=[]
                )
                component_db_model.save()

            components.append(component)

        for workflow_component in definition.components:
            component: Component = [c for c in components if c.id == workflow_component.component_id][0]

            if not component:
                raise NoLabsException(ErrorCodes.component_not_found)

            # Check that component mappings exist
            for mapping in workflow_component.mappings:
                source_component: Component = [c for c in components if c.id == mapping.source_component_id][0]

                if source_component.id not in component.previous_component_ids:
                    component.add_previous(component_id=source_component.id)

                component.try_map_property(component=source_component,
                                           path_from=mapping.source_path,
                                           target_path=mapping.target_path)

            for default in workflow_component.defaults:
                component.try_set_default(default.target_path, value=default.value)

        # executor = Executor(self.logger)
        # await executor.execute(components=components)

    def fetch_workflow(self, workflow_id: uuid.UUID) -> WorkflowDbModel:
        return WorkflowDbModel.objects.with_id(workflow_id)

    def fetch_component(self, component_id: uuid.UUID) -> ComponentDbModel:
        return ComponentDbModel.objects.with_id(component_id)


class StartWorkflowComponentFeature:
    available_components: Dict[str, Type[Component]]
    logger: logging.Logger

    def __init__(self, available_components: Dict[str, Type[Component]],
                 logger: logging.Logger = logging.getLogger('nolabs')):
        self.available_components = available_components
        self.logger = logger

    async def handle(self, request: StartWorkflowComponentRequest):
        db_model: WorkflowDbModel = WorkflowDbModel.objects.with_id(request.workflow_id)
        experiment = db_model.experiment

        workflow_schema = db_model.get_workflow_definition()

        if not workflow_schema.valid:
            raise NoLabsException(ErrorCodes.invalid_workflow_definition,
                                  f'Run {UpdateWorkflowSchemaFeature.__name__} to get schema errors')

        components: List[Component] = []

        component: Component
        for workflow_component in workflow_schema.components:
            component_db_model: ComponentDbModel = ComponentDbModel.objects.with_id(
                workflow_component.component_id)
            if component_db_model:
                component: Component = component_db_model.get_component()
            else:
                id = workflow_component.component_id
                component: Component = self.available_components[workflow_component.name](id=id,
                                                                                          experiment=experiment)
                component_db_model = ComponentDbModel.create(
                    id=id,
                    workflow=db_model,
                    component=component,
                    jobs=[]
                )
                component_db_model.save()

            components.append(component)

        for workflow_component in workflow_schema.components:
            component: Component = [c for c in components if c.id == workflow_component.component_id][0]

            if not component:
                raise NoLabsException(ErrorCodes.component_not_found)

            # Check that component mappings exist
            for mapping in workflow_component.mappings:
                source_component: Component = [c for c in components if c.id == mapping.source_component_id][0]

                if source_component.id in component.previous_component_ids:
                    component.add_previous(component_id=source_component.id)

                component.try_map_property(component=source_component,
                                           path_from=mapping.source_path,
                                           target_path=mapping.target_path)

            for default in workflow_component.defaults:
                component.try_set_default(default.target_path, value=default.value)

        # executor = Executor(logger=self.logger)
        # component = [c for c in components if c.id == request.component_id][0]
        # await executor.execute(components=components, specific=component)


class GetComponentStateFeature:
    async def handle(self, request: GetComponentStateRequest) -> GetComponentStateResponse:
        component: ComponentDbModel = ComponentDbModel.objects.with_id(request.component_id)

        if not component:
            return GetComponentStateResponse(
                input_dict={},
                output_dict={},
                job_ids=[],
                jobs_errors=[],
                last_exceptions=[],
                input_property_errors=[]
            )

        job_errors = []

        job: Job
        for job in Job.objects(id__in=[j.id for j in component.jobs]):
            for error in job.input_errors():
                job_errors.append(JobErrorResponse(
                    msg=error.message,
                    job_id=job.id
                ))

        c = component.get_component()

        return GetComponentStateResponse(
            input_dict=c.input_value_dict,
            output_dict=c.output_value_dict,
            job_ids=[j.iid.value for j in component.jobs],
            jobs_errors=job_errors,
            input_property_errors=[
                InputPropertyErrorResponse(
                    loc=e.loc,
                    msg=e.msg
                )
                for e in component.input_property_errors
            ],
            last_exceptions=[]
        )


class ResetWorkflowFeature:
    # todo change
    async def handle(self, request: ResetWorkflowRequest):
        workflow: WorkflowDbModel = WorkflowDbModel.objects.with_id(request.workflow_id)

        for component in workflow.get_workflow_definition().components:
            db_model: ComponentDbModel = ComponentDbModel.objects.with_id(component.component_id)

            db_model.set_component({})

            db_model.input_parameter_dict = {}
            db_model.output_parameter_dict = {}
            db_model.input_property_errors = []

            db_model.save()
