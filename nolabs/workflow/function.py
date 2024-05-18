__all__ = [
    'PropertyValidationError',
    'PythonFunction',
    'Properties'
]

import asyncio
import inspect
import uuid
from abc import ABC, abstractmethod
from asyncio import Task
from typing import Optional, List, Any, Type, Dict, Union, get_type_hints, TypeVar, Generic

from pydantic import ValidationError, BaseModel

from nolabs.workflow.component import Component
from nolabs.workflow.exceptions import WorkflowException
from nolabs.workflow.properties import ParameterSchema, Property, PropertyValidationError


class Function(ABC):
    id: str
    name: str
    signature: str

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

    @abstractmethod
    def try_map_property(self, function: 'Function', path_from: List[str], path_to: List[str]) -> Optional[
        PropertyValidationError]:
        ...

    @abstractmethod
    def set_properties_from_previous(self):
        ...

    @property
    @abstractmethod
    def unmapped_properties(self) -> List[PropertyValidationError]:
        ...

    @abstractmethod
    def validate_output(self) -> List[PropertyValidationError]:
        ...


TInput = TypeVar('TInput', bound=BaseModel)
TOutput = TypeVar('TOutput', bound=BaseModel)
TParameter = TypeVar('TParameter', bound=BaseModel)


class PythonParameter(Generic[TParameter]):
    _type: Type[TParameter]
    _dictionary: Dict[str, Any] = {}
    schema: ParameterSchema

    def __init__(self, type: Type[TParameter], dictionary: Union[Dict[str, Any], None] = None):
        if dictionary is None:
            dictionary = {}

        if not issubclass(type, BaseModel):
            raise WorkflowException(f"Parameter of type {type} is not a pydantic @dataclass")

        if dictionary and not self.validate_dictionary(dictionary):
            raise WorkflowException(f"Input value type {type} is not compatible with type")

        self._dictionary = dictionary
        self._type = type
        self.schema = ParameterSchema.get_instance(cls=type)

    def validate(self) -> List[PropertyValidationError]:
        return self.validate_dictionary(self._dictionary)

    def validate_dictionary(self, dictionary: Dict[str, Any]) -> List[PropertyValidationError]:
        try:
            _ = self._type(**dictionary)
        except ValidationError as e:
            return [
                PropertyValidationError(
                    msg=error['msg'],
                    loc=error['loc']  # type: ignore
                )
                for error in e.errors()
            ]

        return []

    def set_instance(self, instance: TParameter):
        if not self.is_pydantic_type(instance):
            raise WorkflowException(f'Cannot set instance other than {str(self._type)}')

        value = instance.dict()

        errors = self.validate_dictionary(dictionary=value)
        if errors:
            raise WorkflowException(', '.join([str(err) for err in errors]))

        self._dictionary = value

    def set_parameter(self, path: List[str], param: Any):
        current_level = self._dictionary
        for key in path[:-1]:
            if key not in current_level:
                current_level[key] = {}
            current_level = current_level[key]
        current_level[path[-1]] = param

    def get_parameter(self, path: List[str]) -> Any:
        current_level = self._dictionary
        for key in path[:-1]:
            if key not in current_level:
                current_level[key] = {}
            current_level = current_level[key]
        return current_level[path[-1]]

    @property
    def instance(self) -> TParameter:
        errors = self.validate_dictionary(self._dictionary)

        if errors:
            raise WorkflowException(', '.join([str(err) for err in errors]))

        return self._type(**self._dictionary)

    def reset_value(self):
        self._dictionary = {}

    @property
    def dict(self) -> Dict[str, Any]:
        return self._dictionary

    @staticmethod
    def is_pydantic_type(t: Any) -> bool:
        return issubclass(type(t), BaseModel) or '__pydantic_post_init__' in t.__dict__


class PythonFunction(Function, Generic[TInput, TOutput]):
    _execution_timeout: int = 3600
    _component: Component[TInput, TOutput]

    _input: PythonParameter[TInput]
    _output: PythonParameter[TOutput]

    title: str

    last_exception: Optional[Exception] = None

    previous: List['Function'] = []

    def __init__(self,
                 component: Component[TInput, TOutput],
                 execution_timeout: int = 3600):
        if not isinstance(component, Component):
            raise WorkflowException(f'Component must be instance of {Component}')

        if not inspect.iscoroutinefunction(component.start):
            raise WorkflowException(f'Function {str(component.start)} must be a coroutine')

        self._component = component
        self._execution_timeout = execution_timeout

        name = self._parse_function_name()
        signature = f'def {name}{str(inspect.signature(component.start))}'

        self.id = str(uuid.uuid4())
        self.signature = signature
        self._input = self._parse_input_parameter()
        self._output = self._parse_return_parameter()
        self.title = component.title

    async def execute(self):
        input_validation_errors = self._input.validate()
        if input_validation_errors:
            raise WorkflowException(', '.join([str(err) for err in input_validation_errors]))

        try:
            self._output.reset_value()

            result = await asyncio.create_task(
                asyncio.wait_for(self._component.start(self._input.instance), self._execution_timeout))

            self._output.set_instance(result)
        except Exception as e:
            if not isinstance(e, asyncio.CancelledError):
                self.last_exception = e

    async def terminate(self, timeout: int = 10):
        await asyncio.wait_for(self._component.stop(), timeout=timeout)

    def _parse_input_parameter(self) -> PythonParameter:
        sig = inspect.signature(self._component.start)

        if len(sig.parameters) != 1:
            raise WorkflowException('Your function must have exactly one pydantic @dataclass input parameter')

        parameter = list(sig.parameters.values())[0]

        hints = get_type_hints(self._component.start)
        param_type = hints[parameter.name]

        if not PythonParameter.is_pydantic_type(param_type):
            raise WorkflowException(f'Input parameter with type {param_type} is not a pydantic {str(BaseModel)}')

        return PythonParameter(
            type=param_type
        )

    def set_input(self, instance: TInput):
        self._input.set_instance(instance=instance)

    def get_output(self) -> TOutput:
        return self._output.instance

    def validate_input(self):
        return self._input.validate()

    def _parse_function_name(self) -> str:
        return type(self._component).__name__

    def _parse_return_parameter(self) -> PythonParameter:
        sig = inspect.signature(self._component.start)

        output_param = sig.return_annotation

        if not PythonParameter.is_pydantic_type(output_param):
            raise WorkflowException(
                f'Input parameter with type {output_param} is not a pydantic {str(BaseModel)}')

        return PythonParameter(
            type=output_param
        )

    def remove_previous(self, function: 'PythonFunction'):
        if function not in self.previous:
            return

        self.previous.remove(function)

        if not self._input.schema.properties.items():
            return

        for _, property in self._input.schema.properties.items():
            if property.source_function_id == function.id:
                property.unmap()

    def add_previous(self, function: 'PythonFunction'):
        if function in self.previous:
            return

        self.previous.append(function)

    def try_map_property(self, function: 'Function', path_from: List[str], path_to: List[str]) -> Optional[
        PropertyValidationError]:
        if function not in self.previous:
            raise WorkflowException(f'Cannot map parameter {path_to} for unmapped function {function.id}')

        if not isinstance(function, PythonFunction):
            raise ValueError(f'Function is not a {PythonFunction}')  # TODO change later

        return self._input.schema.try_set_mapping(
            source_schema=function._output.schema,
            function_id=function.id,
            path_from=path_from,
            path_to=path_to
        )

    def set_properties_from_previous(self):
        for function in self.previous:
            if not isinstance(function, PythonFunction):
                raise ValueError(f'Function is not a {PythonFunction}')  # TODO change later
            for prop in self._input.schema.mapped_properties:
                if prop.source_function_id == function.id:
                    input_parameter = function._output.get_parameter(prop.path_from)
                    self._input.set_parameter(path=prop.path_to, param=input_parameter)

    @property
    def input_properties(self) -> Dict[str, Property]:
        if not self._input.schema.properties:
            return {}

        return self._input.schema.properties

    @property
    def output_properties(self) -> Dict[str, Property]:
        if not self._output.schema.properties:
            return {}

        return self._output.schema.properties

    @property
    def unmapped_properties(self) -> List[PropertyValidationError]:
        result = []
        for prop in self._input.schema.unmapped_properties:
            result.append(PropertyValidationError(
                msg='Unmapped property',
                loc=[prop.title]  # type: ignore
            ))
        return result

    def validate_output(self) -> List[PropertyValidationError]:
        return self._output.validate()

    @property
    def component_id(self) -> str:
        return self._component.id

    def load
