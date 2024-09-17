__all__ = [
    "GetJobStateFeature",
    "GetComponentStateFeature",
    "StartWorkflowFeature",
    "ResetWorkflowFeature",
    "AllWorkflowSchemasFeature",
    "GetWorkflowSchemaFeature",
    "StartWorkflowComponentFeature",
    "UpdateWorkflowSchemaFeature",
    "CreateWorkflowSchemaFeature",
    "DeleteWorkflowSchemaFeature",
]

import uuid
from typing import List, Optional
from uuid import UUID

from domain.exceptions import ErrorCodes, NoLabsException
from infrastructure.log import logger
from markdown.extensions import extra
from prefect import get_client
from prefect.client.schemas import StateType
from prefect.client.schemas.objects import TERMINAL_STATES

from nolabs.application.workflow.api.api_models import (
    AllWorkflowSchemasResponse, ComponentStateEnum, GetComponentRequest,
    GetComponentResponse, GetJobRequest, GetJobState, JobStateEnum,
    PropertyErrorResponse, ResetWorkflowRequest, StartWorkflowComponentRequest)
from nolabs.application.workflow.api.mappings import map_property
from nolabs.application.workflow.api.schema import (ComponentSchema,
                                                    ComponentSchemaTemplate,
                                                    WorkflowSchema)
from nolabs.application.workflow.component import (Component,
                                                   ComponentTypeFactory)
from nolabs.application.workflow.dag import PrefectDagExecutor
from nolabs.application.workflow.data import (ComponentData,
                                              ExperimentWorkflowRelation,
                                              JobRunData, WorkflowData)


class DeleteWorkflowSchemaFeature:
    async def handle(self, id: UUID):
        extra = {"workflow_id": id}

        logger.info("Deleting workflow", extra=extra)

        try:
            WorkflowData.objects.with_id(id).delete()

            logger.info("Workflow deleted successfully", extra=extra)
        except Exception as e:
            if e is NoLabsException:
                raise e
            raise NoLabsException(ErrorCodes.delete_workflow_failed) from e


class AllWorkflowSchemasFeature:
    async def handle(self, experiment_id: UUID) -> AllWorkflowSchemasResponse:

        extra = {"experiment_id": experiment_id}

        logger.info("Get all workflow schemas", extra=extra)

        try:
            relation: ExperimentWorkflowRelation = ExperimentWorkflowRelation.objects(
                experiment=experiment_id
            ).first()

            response = AllWorkflowSchemasResponse(
                ids=[relation.workflow.id] if relation else []
            )

            logger.info("Get all workflow schemas success", extra=extra)

            return response
        except Exception as e:
            if e is NoLabsException:
                raise e
            raise NoLabsException(ErrorCodes.get_all_workflows_failed) from e


class CreateWorkflowSchemaFeature:
    async def handle(self, experiment_id: UUID) -> WorkflowSchema:
        extra = {"experiment_id": experiment_id}

        logger.info("Staring workflow creation", extra=extra)

        try:
            id = uuid.uuid4()

            component_templates: List[ComponentSchemaTemplate] = []
            component_schemas: List[ComponentSchema] = []

            for name, component_type in ComponentTypeFactory.enumerate():
                component = component_type(id=uuid.uuid4())

                output_schema = component.output_schema
                output_parameters = {
                    name: map_property(prop, output_schema)
                    for name, prop in output_schema.properties.items()
                }

                input_schema = component.input_schema
                input_parameters = {
                    name: map_property(prop, input_schema)
                    for name, prop in input_schema.properties.items()
                }

                component_templates.append(
                    ComponentSchemaTemplate(
                        name=component.name,
                        input=input_parameters,
                        output=output_parameters,
                        description=component.description,
                    )
                )

            schema = WorkflowSchema(
                workflow_id=id,
                error=None,
                component_templates=component_templates,
                components=component_schemas,
            )

            state = WorkflowData.create(id=id, schema=schema.dict())
            state.save(cascade=True)

            ExperimentWorkflowRelation.create(
                experiment_id=experiment_id, workflow_id=id
            ).save()

            logger.info("Created workflow with id %s", id, extra=extra)

            return schema
        except Exception as e:
            if e is NoLabsException:
                raise e
            raise NoLabsException(ErrorCodes.create_workflow_failed) from e


class GetWorkflowSchemaFeature:
    async def handle(self, id: UUID) -> Optional[WorkflowSchema]:
        extra = {"workflow_id": id}

        logger.info("Get workflow schema", extra=extra)

        try:
            data: WorkflowData = WorkflowData.objects.with_id(id)

            if not data:
                return None

            logger.info("Get workflow schema success", extra=extra)

            return WorkflowSchema(**data.schema)
        except Exception as e:
            if e is NoLabsException:
                raise e
            raise NoLabsException(ErrorCodes.get_workflow_schema_failed) from e


class UpdateWorkflowSchemaFeature:
    async def handle(self, schema: WorkflowSchema) -> WorkflowSchema:
        extra = {
            "workflow_id": schema.workflow_id,
            "components_count": len(schema.components),
        }

        logger.info("Update workflow schema", extra=extra)

        try:
            data = WorkflowData.objects.get(id=schema.workflow_id)

            components: List[Component] = []

            for wf_component in schema.components:
                component = ComponentTypeFactory.get_type(wf_component.name)(
                    id=wf_component.component_id
                )

                components.append(component)

            # Validate graph
            schema_cyclic = self._is_cyclic(graph=components)

            if schema_cyclic:
                schema.error = "Workflow must be acyclic"
                schema.valid = False
                data.schema = schema.dict()
                data.save(cascade=True)

                logger.info(
                    "Update workflow schema success, schema is cyclic", extra=extra
                )
                return schema

            schema_valid = True

            for workflow_component in schema.components:
                component = [
                    c for c in components if c.id == workflow_component.component_id
                ][0]

                # Validate mappings
                for mapping in workflow_component.mappings:
                    source_components = [
                        c
                        for c in schema.components
                        if c.component_id == mapping.source_component_id
                    ]
                    if not source_components:
                        raise NoLabsException(ErrorCodes.component_not_found)

                    source_component = [
                        c for c in components if c.id == mapping.source_component_id
                    ][0]

                    component.add_previous(component_id=source_component.id)

                    error = component.try_map_property(
                        component=source_component,
                        path_from=mapping.source_path,
                        target_path=mapping.target_path,
                    )

                    if error:
                        mapping.error = error.msg
                        schema_valid = False

                # Validate defaults
                for default in workflow_component.defaults:
                    error = component.try_set_default(
                        target_path=default.target_path, value=default.value
                    )
                    if error:
                        default.error = error.msg
                        schema_valid = False

                component_data = ComponentData.objects.with_id(component.id)

                if not component_data:
                    component_data = ComponentData.create(
                        id=component.id, workflow=schema.workflow_id
                    )

                component.dump(data=component_data)
                component_data.save()

            schema.valid = schema_valid
            data.schema = schema.dict()
            data.save()

            logger.info(
                "Update workflow schema success",
                extra={**extra, **{"schema_valid": schema.valid}},
            )

            return schema
        except Exception as e:
            if e is NoLabsException:
                raise e
            raise NoLabsException(ErrorCodes.update_workflow_schema_failed) from e

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
        extra = {"workflow_id": id}

        try:
            logger.info("Starting workflow schema", extra=extra)

            relation: ExperimentWorkflowRelation = ExperimentWorkflowRelation.objects(
                workflow=id
            ).first()

            data = WorkflowData.objects.with_id(id)

            schema = WorkflowSchema(**data.schema)

            if not schema.valid:
                raise NoLabsException(ErrorCodes.invalid_workflow_schema)

            components: List[Component] = []

            for schema_component in schema.components:
                data = ComponentData.objects.with_id(schema_component.component_id)
                components.append(Component.restore(data=data))

            components: List[Component] = [
                Component.restore(ComponentData.objects.with_id(c.component_id))
                for c in schema.components
            ]
            executor = PrefectDagExecutor()

            extra = {
                **extra,
                **{
                    "component_ids": [c.id for c in components],
                    "experiment_id": relation.experiment.id,
                },
            }

            logger.info("Workflow schema execute", extra=extra)

            await executor.execute(
                workflow_id=id,
                components=components,
                experiment_id=relation.experiment.id,
            )

            logger.info("Workflow schema executed", extra=extra)
        except Exception as e:
            if e is NoLabsException:
                raise e
            raise NoLabsException(ErrorCodes.start_workflow_schema_failed) from e


class StartWorkflowComponentFeature:
    async def handle(self, request: StartWorkflowComponentRequest):
        try:
            extra = {
                "component_id": request.component_id,
                "workflow_id": request.workflow_id,
            }

            logger.info("Starting component", extra=extra)

            relation: ExperimentWorkflowRelation = ExperimentWorkflowRelation.objects(
                workflow=request.workflow_id
            ).first()

            data = WorkflowData.objects.with_id(request.workflow_id)

            schema = WorkflowSchema(**data.schema)

            if not schema.valid:
                raise NoLabsException(ErrorCodes.invalid_workflow_schema)

            components: List[Component] = [
                Component.restore(ComponentData.objects.with_id(c.component_id))
                for c in schema.components
            ]
            executor = PrefectDagExecutor()

            extra = {**extra, **{"experiment_id": relation.experiment.id}}

            logger.info("Component execution", extra=extra)

            await executor.execute(
                workflow_id=request.workflow_id,
                components=[c for c in components if c == request.component_id],
                experiment_id=relation.experiment.id,
            )

            logger.info("Component execution finished", extra=extra)
        except Exception as e:
            if e is NoLabsException:
                raise e
            raise NoLabsException(ErrorCodes.start_component_failed) from e


class GetComponentStateFeature:
    async def handle(self, request: GetComponentRequest) -> GetComponentResponse:
        try:
            extra = {"component_id": request.id}

            logger.info("Get component state", extra=extra)

            data: ComponentData = ComponentData.objects.with_id(request.id)

            if not data:
                raise NoLabsException(ErrorCodes.component_not_found)

            if not data.flow_run_id:
                raise NoLabsException(ErrorCodes.flow_run_id_not_found)

            job_ids = [
                j.id for j in JobRunData.objects(component=request.id).only("id")
            ]

            async with get_client() as client:
                flow_run = await client.read_flow_run(data.flow_run_id)

            prefect_state = flow_run.state_type

            state = ComponentStateEnum.RUNNING

            if prefect_state in TERMINAL_STATES:
                state = ComponentStateEnum.COMPLETED

            if prefect_state in [
                StateType.FAILED,
                StateType.CANCELLING,
                StateType.CANCELLED,
                StateType.CRASHED,
            ]:
                state = ComponentStateEnum.FAILED

            response = GetComponentResponse(
                id=data.id,
                input_schema=data.input_schema,
                output_schema=data.output_schema,
                input_value_dict=data.input_value_dict,
                output_value_dict=data.output_value_dict,
                previous_component_ids=data.previous_component_ids,
                input_errors=[
                    PropertyErrorResponse(loc=e.loc, msg=e.msg)
                    for e in data.input_errors
                ],
                output_errors=[
                    PropertyErrorResponse(loc=e.loc, msg=e.msg)
                    for e in data.output_errors
                ],
                state=state,
                state_message=flow_run.state.message,
                job_ids=job_ids,
            )

            logger.info("Get component state success", extra=extra)

            return response
        except Exception as e:
            if e is NoLabsException:
                raise e
            raise NoLabsException(ErrorCodes.get_component_state_failed) from e


class GetJobStateFeature:
    async def handle(self, request: GetJobRequest) -> GetJobState:
        try:
            extra = {"job_id": request.job_id}

            logger.info("Get job state", extra=extra)

            job: JobRunData = JobRunData.objects.with_id(request.job_id)

            if not job:
                raise NoLabsException(ErrorCodes.job_not_found)

            async with get_client() as client:
                task_run = await client.read_task_run(task_run_id=job.task_run_id)

            prefect_state = task_run.state_type

            state = JobStateEnum.RUNNING

            if prefect_state in TERMINAL_STATES:
                state = JobStateEnum.COMPLETED

            if prefect_state in [
                StateType.FAILED,
                StateType.CANCELLING,
                StateType.CANCELLED,
                StateType.CRASHED,
            ]:
                state = JobStateEnum.FAILED

            response = GetJobState(
                id=job.id,
                component_id=job.component_id,
                state=state,
                state_message=task_run.state.message,
            )

            extra = {
                **extra,
                **{"component_id": job.component_id, "job_state": state.name},
            }

            logger.info("Get job state success", extra=extra)

            return response
        except Exception as e:
            if e is NoLabsException:
                raise e
            raise NoLabsException(ErrorCodes.get_job_state_failed) from e


class ResetWorkflowFeature:
    async def handle(self, request: ResetWorkflowRequest):
        pass
