__all__ = [
    'PropertyValidationError',
    'Component',
]

import asyncio
import uuid
from abc import abstractmethod
from typing import Optional, List, Any, Type, Dict, Union, TypeVar, Generic, Tuple, get_args

from pydantic import BaseModel
from pydantic.dataclasses import dataclass

from nolabs.refined.domain.models.common import Job, Experiment
from nolabs.workflow.properties import ParameterSchema, Property, PropertyValidationError


@dataclass
class JobValidationError:
    job_id: uuid.UUID
    msg: str


TInput = TypeVar('TInput', bound=BaseModel)
TOutput = TypeVar('TOutput', bound=BaseModel)


def is_pydantic_type(t: Any) -> bool:
    return issubclass(type(t), BaseModel) or '__pydantic_post_init__' in t.__dict__


class Component(Generic[TInput, TOutput]):
    id: uuid.UUID
    name: str
    description: str
    execution_timeout: int
    experiment: Experiment

    _input_schema: ParameterSchema[TInput]
    _output_schema: ParameterSchema[TOutput]

    input_parameter_dict: Dict[str, Any]
    output_parameter_dict: Dict[str, Any]

    _previous: List['Component']

    jobs: List[Job] = []

    def __init__(self, id: uuid.UUID,
                 experiment: Experiment,
                 jobs: Optional[List[Job]] = None,
                 input_parameter_dict: Optional[Dict[str, Any]] = None,
                 output_parameter_dict: Optional[Dict[str, Any]] = None,
                 execution_timeout: int = 3600):
        self.id = id
        self.execution_timeout = execution_timeout

        self.experiment = experiment

        self.input_parameter_dict = {}
        self.output_parameter_dict = {}

        if not jobs:
            self.jobs = []
        else:
            self.jobs = jobs

        print("INPUT PARAM:", input_parameter_dict)

        if not input_parameter_dict:
            self.input_parameter_dict = {}
        else:
            self.input_parameter_dict = input_parameter_dict

        if not output_parameter_dict:
            self.output_parameter_dict = {}
        else:
            self.output_parameter_dict = output_parameter_dict

        self._input_schema = ParameterSchema.get_instance(cls=self._input_parameter_type)
        self._output_schema = ParameterSchema.get_instance(cls=self._output_parameter_type)

        self._previous = []

    async def executing(self) -> bool:
        return False

    async def terminate(self, timeout: int = 10):
        await asyncio.wait_for(self.stop(), timeout=timeout)

    @property
    def output(self) -> TOutput:
        return self._output_parameter_type(**self.output_parameter_dict)

    @output.setter
    def output(self, output_parameter: Union[TOutput, Dict[str, Any]]):
        if isinstance(output_parameter, BaseModel):
            self.output_parameter_dict = output_parameter.dict()
            return

        self.output_parameter_dict = output_parameter

    @property
    def input(self) -> TInput:
        return self._input_parameter_type(**self.input_parameter_dict)

    @property
    def input_dict(self) -> Dict[str, Any]:
        return self.input_parameter_dict

    @property
    def output_dict(self) -> Dict[str, Any]:
        return self.output_parameter_dict

    def validate_input(self):
        return self._input_schema.validate_dictionary(t=self._input_parameter_type,
                                                      dictionary=self.input_parameter_dict)

    def add_previous(self, component: Union['Component', List['Component']]):
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

        if not isinstance(component, Component):
            raise ValueError(f'Component is not a {Component}')  # TODO change later

        return self._input_schema.try_set_mapping(
            source_schema=component._output_schema,
            component_id=component.id,
            path_from=path_from,
            path_to=path_to
        )

    def try_set_default(self, path_to: List[str], value: Any) -> Optional[PropertyValidationError]:
        return self._input_schema.try_set_default(path_to=path_to, value=value, input_type=self._input_parameter_type)

    def set_input_from_previous(self) -> bool:
        """
        returns: True if input was changed
        """

        changed = False

        for prop in self._input_schema.mapped_properties:
            if prop.default:
                path = prop.path_to

                if not prop.path_to:
                    continue

                current_level = self.input_parameter_dict
                for key in path[:-1]:
                    if key not in current_level:
                        current_level[key] = {}
                    current_level = current_level[key]

                if current_level.get(path[-1]) != prop.default:
                    changed = True

                current_level[path[-1]] = prop.default

        for component in self.previous:
            if not isinstance(component, Component):
                raise ValueError(f'Component is not a {Component}')  # TODO change later
            for prop in self._input_schema.mapped_properties:
                if prop.source_component_id == component.id:


                    current_level = component.output_parameter_dict

                    # Find output parameter from output of previous component

                    path = prop.path_from
                    for key in path[:-1]:
                        if key not in current_level:
                            current_level[key] = {}
                        current_level = current_level[key]
                    input_parameter = current_level[path[-1]]

                    # Find and set input parameter for self function

                    path = prop.path_to

                    current_level = self.input_parameter_dict
                    for key in path[:-1]:
                        if key not in current_level:
                            current_level[key] = {}
                        current_level = current_level[key]

                    if path[-1] not in current_level or current_level[path[-1]] != input_parameter:
                        changed = True

                    current_level[path[-1]] = input_parameter

                    continue

        return changed

    @property
    def output_schema(self) -> ParameterSchema:
        return self._output_schema

    @property
    def input_schema(self) -> ParameterSchema:
        return self._input_schema

    @property
    def previous(self) -> List['Component']:
        return self._previous

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
        return self._output_schema.validate_dictionary(t=self._output_parameter_type,
                                                       dictionary=self.output_parameter_dict)

    def _parse_parameter_types(self) -> Tuple[Type[TInput], Type[TOutput]]:
        args = get_args(self._function.__orig_bases__[0])  # type: ignore
        if not args:
            raise ValueError('Instantiate class with generics specified')
        input_parameter_type, output_parameter_type = args
        return input_parameter_type, output_parameter_type

    @abstractmethod
    async def setup_jobs(self):
        ...

    @abstractmethod
    async def prevalidate_jobs(self) -> List[JobValidationError]:
        ...

    # region Function

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

    @abstractmethod
    async def execute(self):
        pass

    # endregion

    def __hash__(self):
        return self.id.__hash__()

    def __eq__(self, other):
        if not isinstance(other, Component):
            return False

        return self.id == other.id
