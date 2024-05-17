from dataclasses import dataclass
from typing import List, Optional

from nolabs.workflow.exceptions import WorkflowException
from nolabs.workflow.function import PythonFunction


@dataclass
class WorkflowValidationError:
    msg: str


class Workflow:
    functions: List[PythonFunction]

    def __init__(self, functions: List[PythonFunction]):
        error = self.validate(functions)
        if error:
            raise WorkflowException(error.msg)

        self.functions = functions

    async def execute(self, terminate: bool = False):
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

    def validate(self, graph: List[PythonFunction]) -> Optional[WorkflowValidationError]:
        if not graph:
            return WorkflowValidationError(
                msg='Execution graph is empty'
            )

        if self.is_cyclic(graph):
            return WorkflowValidationError(
                msg='Execution graph is cyclic'
            )

        return None

    def is_cyclic(self, graph):
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
