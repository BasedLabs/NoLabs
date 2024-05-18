from dataclasses import dataclass
from typing import List, Optional

from nolabs.workflow.component import Component
from nolabs.workflow.exceptions import WorkflowException
from nolabs.workflow.function import PythonFunction
from nolabs.workflow.workflow_schema import WorkflowSchemaModel


@dataclass
class WorkflowValidationError:
    msg: str


class Workflow:
    running: bool
    functions: List[PythonFunction]

    def __init__(self, functions: List[PythonFunction]):
        error = Workflow.validate_graph(functions)
        if error:
            raise WorkflowException(error.msg)

        self.functions = functions

    async def execute(self, terminate: bool = False):
        try:
            self.running = True

            for function in self.functions:
                if function.unmapped_properties and function.validate_input():
                    raise WorkflowException(f'There were unmapped properties: {function.unmapped_properties}')

            if terminate:
                for function in self.functions:
                    await function.terminate()

            executed: List[PythonFunction] = []

            async def execute(function: PythonFunction):
                nonlocal executed
                for previous_function in function.previous:
                    if previous_function.validate_output():
                        await execute(previous_function)  # type: ignore # TODO Change later

                    # If we get errors anyway - throw exception
                    if previous_function.validate_output():
                        raise WorkflowException(f'Cannot execute function {function.id}')

                function.set_properties_from_previous()

                await function.execute()
                executed.append(function)

            for function in self.functions:
                if function not in executed:
                    await execute(function=function)

        finally:
            self.running = False

    @staticmethod
    def validate_graph(graph: List[PythonFunction]) -> Optional[WorkflowValidationError]:
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
    def _find_function(cls, functions: List[PythonFunction], component_id: str) -> Optional[PythonFunction]:
        for function in functions:
            if function.component_id == component_id:
                return function
        return None

    @classmethod
    def validate_schema(
            cls,
            workflow_schema: WorkflowSchemaModel,
            components: List[Component]
    ) -> bool:
        functions = [
            PythonFunction(component=component) for component in components
        ]

        graph_validation_result = Workflow.validate_graph(graph=functions)
        workflow_schema.error = graph_validation_result.msg

        if Workflow.validate_graph(graph=functions):
            return False

        for workflow_component in workflow_schema.workflow_components:
            function = cls._find_function(component_id=workflow_component.component_id, functions=functions)

            # Check that component exists
            if not function:
                return False

            # Check that function mappings exist
            for mapping in workflow_component.mappings:
                source_function = cls._find_function(component_id=mapping.source_component_id, functions=functions)

                if not source_function:
                    return False

                function.add_previous(function=source_function)
                error = function.try_map_property(function=source_function,
                                                  path_from=mapping.source_path,
                                                  path_to=mapping.target_path)

                if error:
                    return False

        return True

    @staticmethod
    def create_from_schema(workflow_schema: WorkflowSchemaModel, components: List[Component]) -> 'Workflow':
        if not workflow_schema.valid:

    @classmethod
    def set_schema_errors(
            cls,
            workflow_schema: WorkflowSchemaModel,
            components: List[Component]):
        component_not_found_template = lambda c_id: f'Component with id {workflow_component.component_id} not found'

        schema_valid = True
        functions = [
            PythonFunction(component=component) for component in components
        ]

        graph_validation_result = Workflow.validate_graph(graph=functions)
        workflow_schema.error = graph_validation_result.msg

        if Workflow.validate_graph(graph=functions):
            workflow_schema.error = 'Workflow must be acyclic'
            schema_valid = False

        for workflow_component in workflow_schema.workflow_components:
            function = cls._find_function(component_id=workflow_component.component_id, functions=functions)

            # Check that component exists
            if not function:
                workflow_component.error = component_not_found_template(workflow_component.component_id)
                schema_valid = False
                continue

            # Check that function mappings exist
            for mapping in workflow_component.mappings:
                source_function = cls._find_function(component_id=mapping.source_component_id, functions=functions)

                if not source_function:
                    mapping.error = f'Component with id {mapping.source_component_id} does not exist'
                    schema_valid = False
                    continue

                function.add_previous(function=source_function)

                error = function.try_map_property(function=source_function,
                                                  path_from=mapping.source_path,
                                                  path_to=mapping.target_path)

                if error:
                    mapping.error = error.msg
                    schema_valid = False

        workflow_schema.valid = schema_valid
        return schema_valid

    @staticmethod
    def is_cyclic(graph):
        visited: Set[WorkflowGraphNode] = set()  # type: ignore
        recursion_stack = set()

        def dfs(vertex: PythonFunction):
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
