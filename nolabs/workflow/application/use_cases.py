__all__ = [
    "GetJobStateFeature",
    "GetComponentStateFeature",
    "StartWorkflowFeature",
    "GetWorkflowSchemaFeature",
    "StartWorkflowComponentFeature",
    "UpdateWorkflowSchemaFeature",
    "CreateWorkflowSchemaFeature",
]

import uuid
from typing import List, Optional
from uuid import UUID

from domain.exceptions import NoLabsException, ErrorCodes
from domain.models.common import Job, ComponentData, Experiment
from infrastructure.log import nolabs_logger as logger
from infrastructure.redis_client_factory import redlock
from workflow.application.api_models import (
    ComponentStateEnum,
    GetComponentRequest,
    GetComponentResponse,
    GetJobRequest,
    GetJobState,
    JobStateEnum,
    PropertyErrorResponse,
    ResetWorkflowRequest,
    StartWorkflowComponentRequest,
)
from workflow.core.states import TERMINAL_STATES, ControlStates
from workflow.application.mappings import map_property
from workflow.application.schema import (
    ComponentSchema,
    ComponentSchemaTemplate,
    WorkflowSchema,
)
from workflow.core.component import Component, ComponentTypeFactory
from workflow.core.graph import GraphExecutionNode


class CreateWorkflowSchemaFeature:
    async def handle(self, experiment_id: UUID) -> WorkflowSchema:
        extra = {"experiment_id": experiment_id}

        logger.info("Staring workflow creation", extra=extra)

        try:
            id = uuid.uuid4()

            component_templates: List[ComponentSchemaTemplate] = []
            component_schemas: List[ComponentSchema] = []

            experiment = Experiment.objects.with_id(experiment_id)

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
                experiment_id=id,
                error=None,
                component_templates=component_templates,
                components=component_schemas,
            )

            experiment.set_schema(schema=schema.model_dump())
            experiment.save(cascade=True)

            logger.info("Saved schema to experiment %s", id, extra=extra)

            return schema
        except Exception as e:
            if isinstance(e, NoLabsException):
                raise e
            raise NoLabsException(ErrorCodes.create_workflow_failed) from e


class GetWorkflowSchemaFeature:
    async def handle(self, experiment_id: UUID) -> Optional[WorkflowSchema]:
        extra = {"experiment_id": experiment_id}

        logger.info("Get workflow schema", extra=extra)

        try:
            data: Experiment = Experiment.objects.with_id(experiment_id)

            if not data or not data.schema:
                return None

            logger.info("Get workflow schema success", extra=extra)

            return WorkflowSchema(**data.schema)
        except Exception as e:
            if isinstance(e, NoLabsException):
                raise e
            raise NoLabsException(ErrorCodes.get_workflow_schema_failed) from e


class UpdateWorkflowSchemaFeature:
    async def handle(self, schema: WorkflowSchema) -> WorkflowSchema:
        extra = {
            "experiment_id": schema.experiment_id,
            "components_count": len(schema.components),
        }

        logger.info("Update workflow schema", extra=extra)

        try:
            data: Experiment = Experiment.objects.with_id(schema.experiment_id)

            await self._ensure_workflow_not_running(schema.experiment_id)

            components = self._fetch_components(schema=schema)

            # Validate graph
            schema_cyclic = self._is_cyclic(graph=components)

            if schema_cyclic:
                schema.error = "Workflow must be acyclic"
                schema.valid = False
                data.schema = schema.model_dump()
                data.save(cascade=True)

                logger.info(
                    "Update workflow schema success, schema is cyclic",
                    extra=extra
                )
                return schema

            # delete component
            old_comp: ComponentData
            for old_comp in ComponentData.objects(experiment=schema.experiment_id):
                found = False
                for c in components:
                    if c.id == old_comp.id:
                        found = True
                        break
                if not found:
                    await old_comp.delete()

            schema_valid = True

            workflow_component: ComponentSchema
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
                        id=component.id, experiment=schema.experiment_id
                    )

                component.dump(data=component_data)
                component_data.save()

            schema.valid = schema_valid
            data.set_schema(schema=schema.model_dump())
            data.save()

            logger.info(
                "Update workflow schema success",
                extra={**extra, **{"schema_valid": schema.valid}},
            )

            return schema
        except Exception as e:
            if isinstance(e, NoLabsException):
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

    def _fetch_components(self, schema: WorkflowSchema) -> List[Component]:
        components: List[Component] = []

        for wf_component in schema.components:
            component = ComponentTypeFactory.get_type(wf_component.name)(
                id=wf_component.component_id
            )

            components.append(component)

        return components

    async def _ensure_workflow_not_running(self, experiment_id: uuid.UUID):
        graph = GraphExecutionNode(experiment_id=experiment_id)
        return not await graph.started()


class StartWorkflowFeature:
    async def handle(self, experiment_id: UUID):
        extra = {"experiment_id": id}

        try:
            logger.info("Starting workflow schema", extra=extra)

            experiment: Experiment = Experiment.objects.get(id=experiment_id)

            schema = WorkflowSchema(**experiment.schema)

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

            extra = {
                **extra,
                **{
                    "component_ids": [c.id for c in components],
                    "experiment_id": experiment.id,
                },
            }

            logger.info("Workflow schema execute", extra=extra)

            with redlock(key=str(experiment_id)):
                await self._ensure_can_start(experiment_id=experiment_id)
                await self.start_workflow(experiment_id=experiment.id, components_graph=components)
            logger.info("Workflow schema executed", extra=extra)
        except Exception as e:
            if isinstance(e, NoLabsException):
                raise e
            raise NoLabsException(ErrorCodes.start_workflow_failed) from e

    async def _ensure_can_start(self, experiment_id: uuid.UUID):
        graph = GraphExecutionNode(experiment_id=experiment_id)
        if not await graph.started() or not await graph.can_execute():
            raise NoLabsException(ErrorCodes.start_workflow_failed, "Cannot schedule this workflow")

    async def start_workflow(self, experiment_id: uuid.UUID, components_graph: List[Component]):
        graph = GraphExecutionNode(experiment_id=experiment_id)

        if not await graph.can_schedule():
            raise NoLabsException(ErrorCodes.start_workflow_failed, message="Cannot")

        await graph.schedule(components=components_graph)


class StartWorkflowComponentFeature:
    async def handle(self, request: StartWorkflowComponentRequest):
        return
        #try:
        #    extra = {
        #        "component_id": request.component_id,
        #        "experiment_id": request.experiment_id,
        #    }
#
        #    logger.info("Starting component", extra=extra)
#
        #    experiment: Experiment = Experiment.objects.with_id(request.experiment_id)
#
        #    schema = WorkflowSchema(**experiment.schema)
#
        #    if not schema.valid:
        #        raise NoLabsException(ErrorCodes.invalid_workflow_schema)
#
        #    data = ComponentData.objects.with_id(request.component_id)
#
        #    await self._ensure_component_not_running(component_id=data.id)
#
        #    flow_run = await Flow.start_component(experiment_id=request.experiment_id, component_id=request.component_id)
#
        #    extra = {**extra, **{"experiment_id": experiment.id}}
#
        #    logger.info("Component execution", extra=extra)
        #except Exception as e:
        #    if isinstance(e, NoLabsException):
        #        raise e
        #    raise NoLabsException(ErrorCodes.start_component_failed) from e

    #async def _ensure_component_not_running(self, component_id: uuid.UUID):
    #    if not await _ready(id=component_id):
    #        raise NoLabsException(ErrorCodes.component_running)


class GetComponentStateFeature:
    async def handle(self, request: GetComponentRequest) -> GetComponentResponse:
        try:
            extra = {"component_id": request.id}

            logger.info("Get component state", extra=extra)

            data: ComponentData = ComponentData.objects.with_id(request.id)

            if not data:
                raise NoLabsException(ErrorCodes.component_not_found)

            job_ids = [j.id for j in Job.objects(component=request.id).only("id")]

            state, message = await self._get_state(component_id=request.id)

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
                state_message=message,
                job_ids=job_ids,
            )

            logger.info("Get component state success", extra=extra)

            return response
        except Exception as e:
            if isinstance(e, NoLabsException):
                raise e
            raise NoLabsException(ErrorCodes.get_component_state_failed) from e

    async def _get_state(self, experiment_id: uuid.UUID, component_id: uuid.UUID) -> (ComponentStateEnum, str):
        component = GraphExecutionNode(experiment_id=experiment_id).get_component_node(component_id=component_id)
        internal_state = await component.get_state()
        state_message = await component.get_message()

        if internal_state == ControlStates.UNKNOWN:
            return ComponentStateEnum.UNKNOWN, state_message

        state = ComponentStateEnum.RUNNING

        if internal_state in TERMINAL_STATES:
            state = ComponentStateEnum.COMPLETED

        if internal_state == ControlStates.CANCELLED:
            state = ComponentStateEnum.CANCELLED

        if internal_state == ControlStates.FAILURE:
            state = ComponentStateEnum.FAILED

        return state, state_message


class GetJobStateFeature:
    async def handle(self, request: GetJobRequest) -> GetJobState:
        try:
            extra = {"job_id": request.job_id}

            logger.info("Get job state", extra=extra)

            job: Job = Job.objects.with_id(request.job_id)

            if not job:
                raise NoLabsException(ErrorCodes.job_not_found, data={"job_id": job.id})

            state, message = await self._get_state(job.id)

            response = GetJobState(
                id=job.id,
                component_id=job.component.id,
                state=state,
                state_message=message,
            )

            extra = {
                **extra,
                **{"component_id": job.component.id, "job_state": state.name},
            }

            logger.info("Get job state success", extra=extra)

            return response
        except Exception as e:
            if isinstance(e, NoLabsException):
                raise e
            raise NoLabsException(ErrorCodes.get_job_state_failed) from e

    async def _get_state(self, experiment_id: uuid.UUID, component_id: uuid.UUID, job_id: uuid.UUID) -> (JobStateEnum, Optional[str]):
        job = GraphExecutionNode(experiment_id=experiment_id).get_job_node(
            component_id=component_id,
            job_id=job_id
        )
        internal_state = await job.get_state()
        state_message = await job.get_message()

        if internal_state == ControlStates.UNKNOWN:
            return ComponentStateEnum.UNKNOWN, state_message

        state = ComponentStateEnum.RUNNING

        if internal_state in TERMINAL_STATES:
            state = ComponentStateEnum.COMPLETED

        if internal_state == ControlStates.CANCELLED:
            state = ComponentStateEnum.CANCELLED

        if internal_state == ControlStates.FAILURE:
            state = ComponentStateEnum.FAILED

        return state, state_message


class StopWorkflowFeature:
    async def handle(self, request: ResetWorkflowRequest):
        pass
