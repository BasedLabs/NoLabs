__all__ = [
    'GetJobStateFeature',
    'GetComponentStateFeature',
    'StartWorkflowFeature',
    'ResetWorkflowFeature',
    'AllWorkflowSchemasFeature',
    'GetWorkflowSchemaFeature',
    'StartWorkflowComponentFeature',
    'UpdateWorkflowSchemaFeature',
    'CreateWorkflowSchemaFeature',
    'DeleteWorkflowSchemaFeature'
]

import uuid
from typing import List
from typing import Optional
from uuid import UUID

from prefect import get_client
from prefect.client.schemas import StateType
from prefect.client.schemas.objects import TERMINAL_STATES

from nolabs.application.workflow.api.api_models import (GetComponentRequest, GetComponentResponse,
                                                        AllWorkflowSchemasResponse, ResetWorkflowRequest,
                                                        StartWorkflowComponentRequest, PropertyErrorResponse,
                                                        ComponentStateEnum, GetJobRequest, GetJobResponse, JobStateEnum)
from nolabs.application.workflow.api.mappings import map_property
from nolabs.application.workflow.api.schema import WorkflowSchema, ComponentSchemaTemplate, ComponentSchema
from nolabs.application.workflow.component import ComponentTypeFactory, Component
from nolabs.application.workflow.dag import PrefectDagExecutor
from nolabs.application.workflow.data import WorkflowData, ComponentData, JobRunData, ExperimentWorkflowRelation
from domain.exceptions import NoLabsException, ErrorCodes


class DeleteWorkflowSchemaFeature:
    async def handle(self, id: UUID):
        WorkflowData.objects.with_id(id).delete()


class AllWorkflowSchemasFeature:
    async def handle(self, experiment_id: UUID) -> AllWorkflowSchemasResponse:
        relation: ExperimentWorkflowRelation = ExperimentWorkflowRelation.objects(experiment=experiment_id).first()

        return AllWorkflowSchemasResponse(
            ids=[relation.workflow.id] if relation else []
        )


class CreateWorkflowSchemaFeature:
    async def handle(self, experiment_id: UUID) -> WorkflowSchema:
        id = uuid.uuid4()

        component_templates: List[ComponentSchemaTemplate] = []
        component_schemas: List[ComponentSchema] = []

        for name, component_type in ComponentTypeFactory.enumerate():
            component = component_type(id=uuid.uuid4())

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

        state = WorkflowData.create(id=id, schema=schema.dict())
        state.save(cascade=True)

        ExperimentWorkflowRelation.create(
            experiment_id=experiment_id,
            workflow_id=id
        ).save()

        return schema


class GetWorkflowSchemaFeature:
    async def handle(self, id: UUID) -> Optional[WorkflowSchema]:
        data: WorkflowData = WorkflowData.objects.with_id(id)

        if not data:
            return None

        return WorkflowSchema(**data.schema)


class UpdateWorkflowSchemaFeature:
    async def handle(self, schema: WorkflowSchema) -> WorkflowSchema:
        data = WorkflowData.objects.get(id=schema.workflow_id)

        components: List[Component] = []

        for wf_component in schema.components:
            component = ComponentTypeFactory.get_type(wf_component.name)(id=wf_component.component_id)

            components.append(component)

        # Validate graph
        if self._is_cyclic(graph=components):
            schema.error = 'Workflow must be acyclic'
            schema.valid = False
            data.schema = schema.dict()
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

            component_data = ComponentData.objects.with_id(component.id)

            if not component_data:
                component_data = ComponentData.create(id=component.id, workflow=schema.workflow_id)

            component.dump(data=component_data)
            component_data.save()

        schema.valid = schema_valid
        data.schema = schema.dict()
        data.save()

        return schema

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


class StartWorkflowFeature:
    async def handle(self, id: UUID):
        relation: ExperimentWorkflowRelation = ExperimentWorkflowRelation.objects(workflow=id).first()

        data = WorkflowData.objects.with_id(id)

        schema = WorkflowSchema(**data.schema)

        if not schema.valid:
            raise NoLabsException(ErrorCodes.invalid_workflow_schema)

        components: List[Component] = []

        for schema_component in schema.components:
            data = ComponentData.objects.with_id(schema_component.component_id)

            #if data:
            components.append(Component.restore(data=data))
           #else:
           #    component_type = ComponentTypeFactory.get_type(schema_component.name)
           #    component = component_type(id=schema_component.component_id)
           #    data = ComponentData.create(id=schema_component.component_id, workflow=id)
           #    component.dump(data=data)
           #    data.save()
           #    components.append(component)


        components: List[Component] = [Component.restore(ComponentData.objects.with_id(c.component_id)) for c in
                                       schema.components]
        executor = PrefectDagExecutor()
        await executor.execute(workflow_id=id,
                               components=components,
                               experiment_id=relation.experiment.id)


class StartWorkflowComponentFeature:
    async def handle(self, request: StartWorkflowComponentRequest):
        relation: ExperimentWorkflowRelation = ExperimentWorkflowRelation.objects(workflow=request.workflow_id).first()

        data = WorkflowData.objects.with_id(request.workflow_id)

        schema = WorkflowSchema(**data.schema)

        if not schema.valid:
            raise NoLabsException(ErrorCodes.invalid_workflow_schema)

        components: List[Component] = [Component.restore(ComponentData.objects.with_id(c.component_id)) for c in
                                       schema.components]
        executor = PrefectDagExecutor()
        await executor.execute(workflow_id=request.workflow_id,
                               components=[c for c in components if c == request.component_id],
                               experiment_id=relation.experiment.id)


class GetComponentStateFeature:
    async def handle(self, request: GetComponentRequest) -> GetComponentResponse:
        data: ComponentData = ComponentData.objects.with_id(request.id)

        if not data:
            raise NoLabsException(ErrorCodes.component_not_found)

        if not data.flow_run_id:
            raise NoLabsException(ErrorCodes.flow_run_id_not_found)

        job_ids = [j.id for j in JobRunData.objects(component=request.id).only('id')]

        async with get_client() as client:
            flow_run = await client.read_flow_run(data.flow_run_id)

        prefect_state = flow_run.state_type

        state = ComponentStateEnum.RUNNING

        if prefect_state in TERMINAL_STATES:
            state = ComponentStateEnum.COMPLETED

        if prefect_state in [StateType.FAILED, StateType.CANCELLING, StateType.CANCELLED, StateType.CRASHED]:
            state = ComponentStateEnum.FAILED

        return GetComponentResponse(
            id=data.id,
            input_schema=data.input_schema,
            output_schema=data.output_schema,
            input_value_dict=data.input_value_dict,
            output_value_dict=data.output_value_dict,
            previous_component_ids=data.previous_component_ids,
            input_errors=[PropertyErrorResponse(loc=e.loc, msg=e.msg) for e in data.input_errors],
            output_errors=[PropertyErrorResponse(loc=e.loc, msg=e.msg) for e in data.output_errors],
            state=state,
            state_message=flow_run.state.message,
            job_ids=job_ids
        )


class GetJobStateFeature:
    async def handle(self, request: GetJobRequest) -> GetJobResponse:
        job: JobRunData = JobRunData.objects.with_id(request.job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        async with get_client() as client:
            task_run = await client.read_task_run(task_run_id=job.task_run_id)

        prefect_state = task_run.state_type

        state = JobStateEnum.RUNNING

        if prefect_state in TERMINAL_STATES:
            state = JobStateEnum.COMPLETED

        if prefect_state in [StateType.FAILED, StateType.CANCELLING, StateType.CANCELLED,
                             StateType.CRASHED]:
            state = JobStateEnum.FAILED

        return GetJobResponse(
            id=job.id,
            component_id=job.component_id,
            state=state,
            state_message=task_run.state.message,
        )


class ResetWorkflowFeature:
    async def handle(self, request: ResetWorkflowRequest):
        pass
