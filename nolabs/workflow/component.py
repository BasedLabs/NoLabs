__all__ = [
    'PropertyValidationError',
    'PythonComponent',
]

import asyncio
import uuid
from abc import ABC, abstractmethod
from typing import Optional, List, Any, Type, Dict, Union, TypeVar, Generic, Tuple, get_args

from mongoengine import Document, UUIDField, ListField
from pydantic import BaseModel

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.workflow.properties import ParameterSchema, Property, PropertyValidationError


class Component(ABC):
    id: uuid.UUID
    name: str
    job_ids: List[uuid.UUID]

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

class PythonComponentDbModel(Document):
    id: uuid.UUID = UUIDField(primary_key=True)
    job_ids: List[uuid.UUID] = ListField(UUIDField)


def is_pydantic_type(t: Any) -> bool:
    return issubclass(type(t), BaseModel) or '__pydantic_post_init__' in t.__dict__


class PythonComponent(Generic[TInput, TOutput], Component):
    _execution_timeout: int = 3600

    _input_schema: ParameterSchema[TInput]
    _output_schema: ParameterSchema[TOutput]

    _input_parameter_dict: Dict[str, Any]
    _output_parameter_dict: Dict[str, Any]

    last_exception: Optional[Exception] = None

    previous: List['Component']

    def __init__(self,
                 id: uuid.UUID,
                 execution_timeout: int = 3600):
        self._execution_timeout = execution_timeout

        self.id = id
        self._input_schema = ParameterSchema.get_instance(cls=self._input_parameter_type)
        self._output_schema = ParameterSchema.get_instance(cls=self._output_parameter_type)

        self._input_parameter_dict = {}
        self._output_parameter_dict = {}

        self.previous = []

    async def execute(self):
        input_validation_errors = self._input_schema.validate_dictionary(self._input_parameter_dict)
        if input_validation_errors:
            raise NoLabsException(ErrorCodes.component_input_invalid, ', '.join([str(err) for err in input_validation_errors]))

        try:
            await asyncio.create_task(
                asyncio.wait_for(self._execute(), self._execution_timeout))
        except Exception as e:
            if not isinstance(e, asyncio.CancelledError):
                self.last_exception = e

    async def terminate(self, timeout: int = 10):
        await asyncio.wait_for(self.stop(), timeout=timeout)

    def set_input(self, instance: TInput):
        if not is_pydantic_type(instance):
            raise ValueError(f'Type must inherit from {BaseModel}')

        value = instance.dict()

        errors = self._input_schema.validate_dictionary(dictionary=value)
        if errors:
            raise NoLabsException(ErrorCodes.component_input_invalid, ', '.join([str(err) for err in errors]))

        self._input_parameter_dict = instance.dict()

    @property
    def output(self) -> TOutput:
        return self._output_parameter_type(**self._output_parameter_dict)

    @output.setter
    def output(self, output_parameter: Union[TOutput, Dict[str, Any]]):
        if isinstance(output_parameter, BaseModel):
            self._output_parameter_dict = output_parameter.dict()
            return

        self._output_parameter_dict = output_parameter

    @property
    def input(self) -> TInput:
        return self._input_parameter_type(**self._input_parameter_dict)

    @property
    def input_dict(self) -> Dict[str, Any]:
        return self._input_parameter_dict

    @property
    def output_dict(self) -> Dict[str, Any]:
        return self._output_parameter_dict

    def validate_input(self):
        return self._input_schema.validate_dictionary(dictionary=self._input_parameter_dict)

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
            raise ValueError(f'Cannot map parameter {path_to} for unmapped component {component.id}')

        if not isinstance(component, PythonComponent):
            raise ValueError(f'Component is not a {PythonComponent}')  # TODO change later

        return self._input_schema.try_set_mapping(
            source_schema=component._output_schema,
            component_id=component.id,
            path_from=path_from,
            path_to=path_to
        )

    def try_set_default(self, path_to: List[str], value: Any) -> Optional[PropertyValidationError]:
        return self._input_schema.try_set_default(path_to=path_to, value=value)

    def set_properties_from_previous(self):
        for prop in self._input_schema.mapped_properties:
            if prop.default:
                path = prop.path_to

                current_level = self._input_parameter_dict
                for key in path[:-1]:
                    if key not in current_level:
                        current_level[key] = {}
                    current_level = current_level[key]
                current_level[path[-1]] = prop.default

        for component in self.previous:
            if not isinstance(component, PythonComponent):
                raise ValueError(f'Component is not a {PythonComponent}')  # TODO change later
            for prop in self._input_schema.mapped_properties:
                if prop.source_component_id == component.id:

                    path = prop.path_from

                    current_level = component._output_parameter_dict

                    # Find output parameter from output of previous component

                    for key in path[:-1]:
                        if key not in current_level:
                            current_level[key] = {}
                        current_level = current_level[key]
                    input_parameter = current_level[path[-1]]

                    # Find and set input parameter for self function

                    path = prop.path_to

                    current_level = self._input_parameter_dict
                    for key in path[:-1]:
                        if key not in current_level:
                            current_level[key] = {}
                        current_level = current_level[key]
                    current_level[path[-1]] = input_parameter
                    continue


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
        return self._output_schema.validate_dictionary(self._output_parameter_dict)

    def _parse_parameter_types(self) -> Tuple[Type[TInput], Type[TOutput]]:
        args = get_args(self._function.__orig_bases__[0])  # type: ignore
        if not args:
            raise ValueError('Instantiate class with generics specified')
        input_parameter_type, output_parameter_type = args
        return input_parameter_type, output_parameter_type

    def _append_job_id(self, job_id: uuid.UUID):
        db_model = PythonComponentDbModel.objects.with_id(self.id)

        if not db_model:
            db_model = PythonComponentDbModel(
                id=self.id,
                job_ids=[]
            )

        db_model.job_ids.append(job_id)
        db_model.save()

    def get_job_ids(self) -> List[uuid.UUID]:
        db_model = PythonComponentDbModel.objects.with_id(self.id)

        if not db_model:
            return []

        return db_model.job_ids

    # region Function

    @abstractmethod
    def _execute(self):
        ...

    @abstractmethod
    async def restore_output(self):
        """Restores output and job_ids"""
        ...

    @property
    @abstractmethod
    def _input_parameter_type(self) -> Type[TInput]:
        ...

    @property
    @abstractmethod
    def _output_parameter_type(self) -> Type[TOutput]:
        ...

    async def stop(self):
        ...

    # endregion
