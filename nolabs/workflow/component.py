__all__ = [
    'PropertyValidationError',
    'PythonComponent',
]

import asyncio
import inspect
import uuid
from abc import ABC, abstractmethod
from typing import Optional, List, Any, Type, Dict, Union, get_type_hints, TypeVar, Generic, Tuple, get_args

from pydantic import ValidationError, BaseModel

from nolabs.workflow.function import PythonFunction
from nolabs.workflow.exceptions import WorkflowException
from nolabs.workflow.properties import ParameterSchema, Property, PropertyValidationError


class Component(ABC):
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
        if not isinstance(other, Component):
            return False

        return self.id == other.id

    @abstractmethod
    def try_map_property(self, component: 'Component', path_from: List[str], path_to: List[str]) -> Optional[
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


def is_pydantic_type(t: Any) -> bool:
    return issubclass(type(t), BaseModel) or '__pydantic_post_init__' in t.__dict__




class PythonComponent(Generic[TInput, TOutput], Component):
    _execution_timeout: int = 3600
    _function: PythonFunction[TInput, TOutput]

    _input_schema: ParameterSchema[TInput]
    _output_schema: ParameterSchema[TOutput]

    last_exception: Optional[Exception] = None

    previous: List['Component']

    def __init__(self,
                 function: PythonFunction[TInput, TOutput],
                 execution_timeout: int = 3600):
        if not isinstance(function, PythonFunction):
            raise WorkflowException(f'Function must be instance of {PythonFunction}')

        self._function = function
        self._execution_timeout = execution_timeout

        if 'execute' not in dir(function):
            raise ValueError('Function must have execute method')

        name = self._parse_function_name()
        signature = f'def {name}{str(inspect.signature(function.execute))}'

        input_type, output_type = self._parse_parameter_types()

        self.id = str(uuid.uuid4())
        self.signature = signature
        self._input_schema = ParameterSchema.get_instance(cls=input_type)
        self._output_schema = ParameterSchema.get_instance(cls=output_type)

        self.previous = []

    async def execute(self):
        input_validation_errors = self._input_schema.validate_dictionary(self._function.input_parameter_dict)
        if input_validation_errors:
            raise WorkflowException(', '.join([str(err) for err in input_validation_errors]))

        try:
            await asyncio.create_task(
                asyncio.wait_for(self._function.execute(), self._execution_timeout))
        except Exception as e:
            if not isinstance(e, asyncio.CancelledError):
                self.last_exception = e

    async def terminate(self, timeout: int = 10):
        await asyncio.wait_for(self._function.stop(), timeout=timeout)

    def set_input(self, instance: TInput):
        if not is_pydantic_type(instance):
            raise WorkflowException(f'Cannot set instance other than {str(TInput)}')

        value = instance.dict()

        errors = self._input_schema.validate_dictionary(dictionary=value)
        if errors:
            raise WorkflowException(', '.join([str(err) for err in errors]))

        self._function.set_input_parameter(instance)

    def get_output(self) -> TOutput:
        return self._function.output_parameter

    def validate_input(self):
        return self._input_schema.validate_dictionary(dictionary=self._function.input_parameter_dict)

    def _parse_function_name(self) -> str:
        return type(self._function).__name__

    def remove_previous(self, component: 'PythonComponent'):
        if component not in self.previous:
            return

        self.previous.remove(component)

        if not self._input_schema.schema.properties.items():
            return

        for _, property in self._input_schema.properties.items():
            if property.source_component_id == component.id:
                property.unmap()

    def add_previous(self, component: Union['PythonComponent', List['PythonComponent']]):
        if isinstance(component, list):
            for c in component:
                if c not in self.previous:
                    self.previous.append(c)

            return

        if component in self.previous:
            return

        self.previous.append(component)

    def try_map_property(self, component: 'Component', path_from: List[str], path_to: List[str]) -> Optional[
        PropertyValidationError]:
        if component not in self.previous:
            raise WorkflowException(f'Cannot map parameter {path_to} for unmapped component {component.id}')

        if not isinstance(component, PythonComponent):
            raise ValueError(f'Component is not a {PythonComponent}')  # TODO change later

        return self._input_schema.try_set_mapping(
            source_schema=component._output_schema,
            component_id=component.id,
            path_from=path_from,
            path_to=path_to
        )

    def set_properties_from_previous(self):
        for component in self.previous:
            if not isinstance(component, PythonComponent):
                raise ValueError(f'Component is not a {PythonComponent}')  # TODO change later
            for prop in self._input_schema.mapped_properties:
                if prop.source_component_id == component.id:

                    path = prop.path_from

                    current_level = component._function.output_parameter_dict

                    # Find output parameter from output of previous component

                    for key in path[:-1]:
                        if key not in current_level:
                            current_level[key] = {}
                        current_level = current_level[key]
                    input_parameter = current_level[path[-1]]

                    # Find and set input parameter for self function

                    path = prop.path_to

                    current_level = self._function.input_parameter_dict
                    for key in path[:-1]:
                        if key not in current_level:
                            current_level[key] = {}
                        current_level = current_level[key]
                    current_level[path[-1]] = input_parameter

    @property
    def input_properties(self) -> Dict[str, Property]:
        if not self._input_schema.properties:
            return {}

        return self._input_schema.properties

    @property
    def output_properties(self) -> Dict[str, Property]:
        if not self._output_schema.properties:
            return {}

        return self._output_schema.properties

    @property
    def unmapped_properties(self) -> List[PropertyValidationError]:
        result = []
        for prop in self._input_schema.unmapped_properties:
            result.append(PropertyValidationError(
                msg='Unmapped property',
                loc=[prop.title]  # type: ignore
            ))
        return result

    def validate_output(self) -> List[PropertyValidationError]:
        return self._output_schema.validate_dictionary(self._function.output_parameter_dict)

    @property
    def component_id(self) -> str:
        return str(self._function.id)

    def _parse_parameter_types(self) -> Tuple[Type[TInput], Type[TOutput]]:
        args = get_args(self._function.__orig_bases__[0])
        if not args:
            raise ValueError('Instantiate class with generics specified')
        input_parameter_type, output_parameter_type = args
        return input_parameter_type, output_parameter_type

    async def restore_function(self):
        self._function.restore()

