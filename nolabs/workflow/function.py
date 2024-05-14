import asyncio
import inspect
import uuid
from abc import ABC, abstractmethod
from asyncio import Task
from dataclasses import dataclass, is_dataclass
from typing import Optional, List, Callable, Set, Any, Type, Dict, Union, Tuple
from uuid import UUID

from pydantic import ValidationError, BaseModel

from nolabs.workflow.schema import Schema, SchemaValidationIssue


@dataclass
class ParameterValidationResult:
    msg: str
    loc: Tuple[str]

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'{self.msg}: {self.loc}'


class Pipe:
    type: Type
    value: Dict[str, Any] = {}
    schema: Schema

    def __init__(self, type: Type, value: Dict[str, Any] = None):
        if not issubclass(type, BaseModel):
            raise ValueError("Parameter is not a pydantic @dataclass")

        if value and not self.validate_value(value):
            raise ValueError("Input value type is not compatible with type")

        self.value = value
        self.schema = Schema.get_schema(cls=type)

    def validate(self) -> List[ParameterValidationResult]:
        return self.validate_value(self.value)

    def validate_value(self, value: Union[Dict[str, Any], Any]) -> List[ParameterValidationResult]:
        try:
            _ = self.type(**value)
        except ValidationError as e:
            return [
                ParameterValidationResult(
                    msg=error.msg,
                    loc=error.loc
                )
                for error in e.errors()
            ]

        return []

    def set_value(self, value: Any):
        errors = self.validate_value(value=value)
        if errors:
            raise ValueError(errors)

        self.value = value

    def set_parameter(self, path: List[str], param: Any):
        current_level = self.value
        for key in path[:-1]:
            if key not in current_level:
                current_level[key] = {}
            current_level = current_level[key]
        current_level[path[-1]] = param

    def get_parameter(self, path: List[str]) -> Any:
        current_level = self.value
        for key in path[:-1]:
            if key not in current_level:
                current_level[key] = {}
            current_level = current_level[key]
        return current_level[path[-1]]

    @property
    def instance(self) -> Type:
        errors = self.validate_value(self.value)

        if errors:
            raise ValueError(errors)

        return self.type(**self.value)

    def reset_value(self):
        self.value = {}

    def to_dict(self) -> Dict[str, Any]:
        return self.value


class Function(ABC):
    id: UUID
    name: str
    signature: str

    input: Pipe
    output: Pipe

    executing: bool = False
    last_exception: Optional[Exception] = None

    previous: List['Function']
    next: List['Function']

    def __init__(self,
                 id: UUID,
                 name: str,
                 signature: str,
                 input: Pipe,
                 output: Pipe):
        self.id = id
        self.name = name
        self.signature = signature
        self.input = input
        self.output = output

    @abstractmethod
    async def execute(self):
        ...

    @abstractmethod
    def terminate(self):
        ...

    def set_previous(self, functions: List['Function']):
        self.previous = functions

    def set_next(self, functions: List['Function']):
        self.next = functions

    def map_parameter(self, function: 'Function', path_from: List[str], path_to: List[str]) -> Optional[
        SchemaValidationIssue]:
        if function not in self.previous:
            raise ValueError('Cannot map parameter for unmapped function')

        return self.input.schema.try_set_mapping(
            source_schema=function.output.schema,
            function_id=function.id,
            path_from=path_from,
            path_to=path_to
        )

    def set_parameters(self, from_function: 'Function'):
        if from_function not in self.previous:
            raise ValueError('Cannot set parameter for unmapped function')

        for prop in self.input.schema.mapped_properties:
            input_parameter = from_function.output.get_parameter(prop.path_from)
            self.input.set_parameter(path=prop.path_to, param=input_parameter)


class PythonFunction(Function):
    _pointer: Callable
    _task: Optional[Task] = None
    _execution_timeout: int = 3600

    def __init__(self,
                 pointer: Callable,
                 execution_timeout: int = 3600):
        if not inspect.iscoroutinefunction(pointer):
            raise ValueError('Function must be a coroutine')

        self._pointer = pointer
        self._execution_timeout = execution_timeout

        name = self._parse_function_name()
        signature = f'def {name}{str(inspect.signature(pointer))}'

        super().__init__(
            id=uuid.uuid4(),
            name=name,
            signature=signature,
            input=self._parse_input_parameter(),
            output=self._parse_return_parameter()
        )

    async def execute(self):
        input_validation_errors = self.input.validate()
        if input_validation_errors:
            raise ValueError(input_validation_errors)

        self.executing = True

        try:
            self.output.reset_value()

            self._task = asyncio.create_task(
                asyncio.wait_for(self._pointer(self.input.instance), self._execution_timeout))

            result = await self._task

            self.output.set_value(result)
        except Exception as e:
            if not isinstance(e, asyncio.CancelledError):
                self.last_exception = e
        finally:
            self._task = None

        self.executing = False

    async def terminate(self, timeout: int = 10):
        if not self._task:
            return

        self._task.cancel()

        async def condition():
            while self._task is not None:
                await asyncio.sleep(0.1)

        await asyncio.wait_for(condition(), timeout=timeout)

    def _parse_input_parameter(self) -> Pipe:
        sig = inspect.signature(self._pointer)

        if len(sig.parameters) != 1:
            raise ValueError('Your function must have exactly one pydantic @dataclass input parameter')

        parameter = list(sig.parameters.values())[0]

        if not is_dataclass(parameter.annotation):
            raise ValueError('Input parameter is not a pydantic @dataclass')

        return Pipe(
            type=parameter.annotation
        )

    def _parse_function_name(self) -> str:
        return self._pointer.__name__

    def _parse_return_parameter(self) -> Pipe:
        sig = inspect.signature(self._pointer)

        return Pipe(
            type=sig.return_annotation
        )

    def __hash__(self):
        return hash(self._pointer)

    def __eq__(self, other) -> bool:
        if not isinstance(other, PythonFunction):
            return False

        return self._pointer == other._pointer


@dataclass
class WorkflowGraphValidationResult:
    valid: bool
    problems: List[str]


class WorkflowContext:
    graph: List[Function]

    def __init__(self, graph: List[Function]):
        validation_result = self.validate(graph)
        if not validation_result.valid:
            raise ValueError('Graph validation issues: ', ', '.join(validation_result.problems))

        self.graph = graph

    async def execute_single(self,
                             function: Function,
                             terminate: bool = False):
        found = False

        def validate_node_in_graph(validation_function: Function):
            global found
            if validation_function == function:
                found = True  # type: ignore

            found = False  # type: ignore

        self._traverse(self.graph[0], validate_node_in_graph, visited=set())

        if not found:
            raise ValueError('Node is not in graph')

        execution_queue: List[Task]

        for previous_node in node.previous:
            if not previous_node.function.output.validate():
                return

        if node.function.executing and terminate:
            await node.function.terminate()

        await node.function.execute()

    async def execute(self, terminate: bool = False):
        if terminate:
            to_terminate: Set[WorkflowGraphNode] = set()

            def _(node: WorkflowGraphNode):
                if node.function.executing and node not in to_terminate:
                    to_terminate.add(node)

            self._traverse(self.graph[0], _, set())

            for node in to_terminate:
                await node.function.terminate()

        head: List[WorkflowGraphNode] = []

        def append_to_queue(node: WorkflowGraphNode):
            if not node.next:
                head.append(node)

        self._traverse(self.graph[0], append_to_queue, set())

        problems = []

        async def execute(node: WorkflowGraphNode):
            """Backwards propagation of execution (from head to tail - dependencies)"""
            nonlocal problems
            input_validation_problems = node.function.input.validate()
            if not input_validation_problems:
                matcher = node.parameters_matcher

                for entry in matcher.entries:
                    for prev_node in node.previous:
                        if prev_node.function == entry.function:

                await node.function.execute()
            else:
                for prev_node in node.previous:
                    await execute(prev_node)
                    errors = prev_node.function.output.validate()
                    if errors:
                        problems += errors
                        return

                for input_parameter in node.function.inputs:
                    for previous_node in node.previous:
                        if input_parameter.type_is_compatible(previous_node.function.returns.type):
                            input_parameter.set_value(previous_node.function.returns.value)

                await node.function.execute()

        for head_node in head:
            await execute(head_node)

    def validate(self, graph: List[Function]) -> WorkflowGraphValidationResult:
        if not graph:
            return WorkflowGraphValidationResult(
                valid=False,
                problems=['Execution graph is empty']
            )

        if self.is_cyclic(graph):
            return WorkflowGraphValidationResult(
                valid=False,
                problems=['Execution graph is cyclic']
            )

        problems: List[str] = []

        def validate_call_chain(function: Function):
            prev_return_params = [function.output for node in node.previous]

            for input_param in function.inputs:
                validated = False
                for ret_param in prev_return_params:
                    if input_param.type_is_compatible(ret_param.type):
                        validated = True

                if not validated:
                    problems.append(
                        f'Input parameter {input_param.name} is not filled for function {node.function.name}')

        self._traverse(node=graph[0], action=validate_call_chain, visited=set())

        return WorkflowGraphValidationResult(
            valid=not problems,
            problems=problems
        )

    def is_cyclic(self, graph):
        visited: Set[WorkflowGraphNode] = set()  # type: ignore
        recursion_stack = set()

        def dfs(vertex: Function):
            if vertex in recursion_stack:
                return True
            if vertex in visited:
                return False

            visited.add(vertex)
            recursion_stack.add(vertex)

            for neighbor in vertex.previous + vertex.next:
                if dfs(neighbor):
                    return True

            recursion_stack.remove(vertex)
            return False

        for vertex in graph:
            if vertex not in visited:
                if dfs(vertex):
                    return True

        return False

    def _traverse(self,
                  node: Function,
                  action: Callable[[Function], None],
                  visited: Set[Function]):
        if node in visited:
            return

        action(node)

        visited.add(node)

        for previous_node in node.previous:
            self._traverse(node=previous_node, action=action, visited=visited)

        for next_node in node.next:
            self._traverse(node=next_node, action=action, visited=visited)
