import asyncio
import inspect
from abc import ABC, abstractmethod
from asyncio import Task
from dataclasses import dataclass
from typing import Optional, List, Callable, Set, Any, Type, Dict, Union, Tuple, get_type_hints

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

    def __init__(self, type: Type, value: Union[Dict[str, Any], None] = None):
        if value is None:
            value = {}

        if not issubclass(type, BaseModel):
            raise ValueError("Parameter is not a pydantic @dataclass")

        if value and not self.get_value_errors(value):
            raise ValueError("Input value type is not compatible with type")

        self.value = value
        self.type = type
        self.schema = Schema.get_schema(cls=type)

    def get_errors(self) -> List[ParameterValidationResult]:
        return self.get_value_errors(self.value)

    def get_value_errors(self, value: Union[Dict[str, Any], Any]) -> List[ParameterValidationResult]:
        try:
            _ = self.type(**value)
        except ValidationError as e:
            return [
                ParameterValidationResult(
                    msg=error['msg'],
                    loc=error['loc']  # type: ignore
                )
                for error in e.errors()
            ]

        return []

    def set_value(self, value: Any):
        if self.is_pydantic_type(value):
            value = value.dict()

        if value is self.type:
            self.value = value.dict()
            return

        errors = self.get_value_errors(value=value)
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
        errors = self.get_value_errors(self.value)

        if errors:
            raise ValueError(errors)

        return self.type(**self.value)

    def reset_value(self):
        self.value = {}

    def to_dict(self) -> Dict[str, Any]:
        return self.value

    @staticmethod
    def is_pydantic_type(t: Type) -> bool:
        return issubclass(type(t), BaseModel) or '__pydantic_post_init__' in t.__dict__


class Function(ABC):
    id: str
    name: str
    signature: str

    input: Pipe
    output: Pipe

    executing: bool = False
    last_exception: Optional[Exception] = None

    previous: List['Function'] = []

    def __init__(self,
                 id: str,
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

    def __hash__(self):
        return self.id.__hash__()

    def __eq__(self, other):
        if not isinstance(other, Function):
            return False

        return self.id == other.id

    def set_previous(self, functions: List['Function']):
        self.previous = functions

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

    def set_parameters_from_previous(self):
        for function in self.previous:
            for prop in self.input.schema.mapped_properties:
                if prop.mapping_function_id == function.id:
                    input_parameter = function.output.get_parameter(prop.path_from)
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
            id=f'{pointer.__module__}.{pointer.__qualname__}',
            name=name,
            signature=signature,
            input=self._parse_input_parameter(),
            output=self._parse_return_parameter()
        )

    async def execute(self):
        input_validation_errors = self.input.get_errors()
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

        hints = get_type_hints(self._pointer)
        param_type = hints[parameter.name]

        if not Pipe.is_pydantic_type(param_type):
            raise ValueError('Input parameter is not a pydantic @dataclass')

        return Pipe(
            type=param_type
        )

    def _parse_function_name(self) -> str:
        return self._pointer.__name__

    def _parse_return_parameter(self) -> Pipe:
        sig = inspect.signature(self._pointer)

        output_param = sig.return_annotation

        if not Pipe.is_pydantic_type(output_param):
            raise ValueError('You can use only pydantic BaseModel inherited class as return value of your function')

        return Pipe(
            type=output_param
        )


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

    async def execute(self, terminate: bool = False):
        if terminate:
            for function in self.graph:
                await function.terminate()

        executed: List[Function] = []

        async def execute(function: Function):
            nonlocal executed
            """Backwards propagation of execution (from head to tail - dependencies)"""
            for previous_function in function.previous:
                if previous_function.output.get_errors():
                    await execute(previous_function)

                # If we get errors anyway - throw exception
                if previous_function.output.get_errors():
                    raise ValueError(f'Cannot execute function {function.id}')

            function.set_parameters_from_previous()

            await function.execute()
            executed.append(function)

        for function in self.graph:
            if function not in executed:
                await execute(function=function)

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

        return WorkflowGraphValidationResult(
            valid=True,
            problems=[]
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

            for neighbor in vertex.previous:
                if dfs(neighbor):
                    return True

            recursion_stack.remove(vertex)
            return False

        for vertex in graph:
            if vertex not in visited:
                if dfs(vertex):
                    return True

        return False
