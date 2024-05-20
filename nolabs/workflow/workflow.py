from dataclasses import dataclass
from typing import List, Optional, Dict
from uuid import UUID

from nolabs.workflow.component import PythonFunction
from nolabs.workflow.exceptions import WorkflowException
from nolabs.workflow.component import PythonComponent
from nolabs.workflow.workflow_schema import WorkflowSchemaModel


@dataclass
class WorkflowValidationError:
    msg: str


class Workflow:
    running: bool
    components: List[PythonComponent]
    executing_workflows: Dict[UUID, 'Workflow']

    def __init__(self, components: List[PythonComponent]):
        error = Workflow.validate_graph(components)
        if error:
            raise WorkflowException(error.msg)

        self.components = components

    async def execute(self,
                      terminate: bool = False,
                      load_parameters: bool = False):
        try:
            self.running = True

            for component in self.components:
                if component.unmapped_properties and component.validate_input():
                    raise WorkflowException(f'There were unmapped properties: {component.unmapped_properties}')

            if terminate:
                for component in self.components:
                    await component.terminate()

            executed: List[PythonComponent] = []

            async def execute(component: PythonComponent):
                nonlocal executed
                for previous_component in component.previous:
                    if previous_component.validate_output():
                        await execute(previous_component)  # type: ignore # TODO Change later

                    # If we get errors anyway - throw exception
                    if previous_component.validate_output():
                        raise WorkflowException(f'Cannot execute component {component.id}')

                component.set_properties_from_previous()

                await component.execute()
                executed.append(component)

            for component in self.components:
                if component not in executed:
                    await execute(component=component)

        finally:
            self.running = False

    @staticmethod
    def validate_graph(graph: List[PythonComponent]) -> Optional[WorkflowValidationError]:
        if not graph:
            return WorkflowValidationError(
                msg='Workflow graph is empty'
            )

        if Workflow.is_cyclic(graph):
            return WorkflowValidationError(
                msg='Workflow graph must be acyclic'
            )

        return None

    @classmethod
    def _find_component(cls, components: List[PythonComponent], component_id: str) -> Optional[PythonComponent]:
        for component in components:
            if component.component_id == component_id:
                return component
        return None

    @classmethod
    def validate_schema(
            cls,
            workflow_schema: WorkflowSchemaModel,
            functions: List[PythonFunction]
    ) -> bool:
        components = [
            PythonComponent(function=function) for function in functions
        ]

        if Workflow.validate_graph(graph=components):
            return False

        for workflow_component in workflow_schema.workflow_components:
            component = cls._find_component(component_id=workflow_component.component_id, components=components)

            # Check that component exists
            if not component:
                return False

            # Check that component mappings exist
            for mapping in workflow_component.mappings:
                source_component = cls._find_component(component_id=mapping.source_component_id, components=components)

                if not source_component:
                    return False

                component.add_previous(component=source_component)
                error = component.try_map_property(component=source_component,
                                                   path_from=mapping.source_path,
                                                   path_to=mapping.target_path)

                if error:
                    return False

        return True

    @classmethod
    def create_from_schema(cls, workflow_schema_model: WorkflowSchemaModel,
                           functions: List[PythonFunction]) -> 'Workflow':
        error = cls.validate_schema(workflow_schema=workflow_schema_model, functions=functions)

        if error:
            raise WorkflowException(
                msg=f'Schema is nov valid. Run {cls.set_schema_errors} to get schema errors'
            )

        components = [
            PythonComponent(function=function) for function in functions
        ]

        for workflow_component in workflow_schema_model.workflow_components:
            component: PythonComponent | None = cls._find_component(component_id=workflow_component.component_id,
                                                                    components=components)

            if not component:
                raise WorkflowException(
                    msg=f'Component with id {workflow_component.component_id} not found'
                )

            # Check that component mappings exist
            for mapping in workflow_component.mappings:
                source_component: PythonComponent | None = cls._find_component(component_id=mapping.source_component_id,
                                                                               components=components)

                if not source_component:
                    raise WorkflowException(
                        msg=f'Component with id {mapping.source_component_id} not found'
                    )

                component.add_previous(component=source_component)
                component.try_map_property(component=source_component,
                                           path_from=mapping.source_path,
                                           path_to=mapping.target_path)

        return Workflow(
            components=components
        )

    @classmethod
    def set_schema_errors(
            cls,
            workflow_schema: WorkflowSchemaModel,
            functions: List[PythonFunction]):
        component_not_found_template = lambda c_id: f'Component with id {c_id} not found'

        schema_valid = True
        components = [
            PythonComponent(function=function) for function in functions
        ]

        graph_validation_error = Workflow.validate_graph(graph=components)

        if graph_validation_error:
            workflow_schema.error = graph_validation_error.msg

        if Workflow.validate_graph(graph=components):
            workflow_schema.error = 'Workflow must be acyclic'
            schema_valid = False

        for workflow_component in workflow_schema.workflow_components:
            component = cls._find_component(component_id=workflow_component.component_id, components=components)

            # Check that component exists
            if not component:
                workflow_component.error = component_not_found_template(workflow_component.component_id)
                schema_valid = False
                continue

            # Check that component mappings exist
            for mapping in workflow_component.mappings:
                source_component = cls._find_component(component_id=mapping.source_component_id, components=components)

                if not source_component:
                    mapping.error = f'Component with id {mapping.source_component_id} does not exist'
                    schema_valid = False
                    continue

                component.add_previous(component=source_component)

                error = component.try_map_property(component=source_component,
                                                   path_from=mapping.source_path,
                                                   path_to=mapping.target_path)

                if error:
                    mapping.error = error.msg
                    schema_valid = False

        workflow_schema.valid = schema_valid
        return schema_valid

    @staticmethod
    def is_cyclic(graph: List[PythonComponent]) -> bool:
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
