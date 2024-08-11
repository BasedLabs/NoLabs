__all__ = [
    'PropertyValidationError'
]

import uuid
from abc import abstractmethod, ABC
from typing import Optional, List, Any, Type, Dict, Union, TypeVar, Generic

from airflow.models import BaseOperator
from airflow.utils.context import Context
from airflow.utils.decorators import apply_defaults
from pydantic import BaseModel
from pydantic.dataclasses import dataclass

from nolabs.application.workflow.components_repository import WorkflowRepository
from nolabs.application.workflow.properties import ParameterSchema, PropertyValidationError
from nolabs.domain.models.common import JobId


@dataclass
class JobValidationError:
    job_id: uuid.UUID
    msg: str


TInput = TypeVar('TInput', bound=BaseModel)
TOutput = TypeVar('TOutput', bound=BaseModel)


def is_pydantic_type(t: Any) -> bool:
    return issubclass(type(t), BaseModel) or '__pydantic_post_init__' in t.__dict__


class Component(ABC, Generic[TInput, TOutput], BaseModel):
    id: uuid.UUID
    experiment_id: uuid.UUID

    job_ids: List[uuid.UUID]

    # region schemas

    output_schema: ParameterSchema[TOutput]
    output_value_dict: Dict[str, Any]

    input_schema: ParameterSchema[TInput]
    input_value_dict: Dict[str, Any]

    previous_component_ids: List[uuid.UUID] = []

    # endregion

    name: str
    description: str

    def create(self, id: uuid.UUID, experiment_id: uuid.UUID):
        c = self.__class__(id=id, experiment_id=experiment_id)
        c.input_schema = ParameterSchema.get_instance(cls=self.input_parameter_type)
        c.output_schema = ParameterSchema.get_instance(cls=self.output_parameter_type)

    @property
    def output_value(self) -> TOutput:
        return self.output_parameter_type(**self.output_value_dict)

    def output_errors(self) -> List[PropertyValidationError]:
        return self.schema.validate_dictionary(t=self.output_parameter_type, dictionary=self.output_value_dict)

    @abstractmethod
    @property
    def input_parameter_type(self) -> Type[TInput]:
        ...

    @abstractmethod
    @property
    def output_parameter_type(self) -> Type[TOutput]:
        ...

    @property
    def input_value(self) -> TInput:
        return self.input_parameter_type(**self.input_value_dict)

    def input_errors(self):
        return self.input_schema.validate_dictionary(t=self.input_parameter_type, dictionary=self.input_value_dict)

    def try_map_property(self, component: 'Component', path_from: List[str], target_path: List[str]) -> Optional[
        PropertyValidationError]:
        return self.input_schema.try_set_mapping(
            source_schema=component.output_schema,
            component_id=component.id,
            path_from=path_from,
            target_path=target_path
        )

    def try_set_default(self, target_path: List[str], value: Any) -> Optional[PropertyValidationError]:
        return self.input_schema.try_set_default(target_path=target_path, value=value,
                                                 input_type=self.input_parameter_type)

    def add_previous(self, component_id: Union[uuid.UUID, List[uuid.UUID]]):
        if isinstance(component_id, list):
            for c in component_id:
                if c not in self.previous_component_ids:
                    self.previous_component_ids.append(c)

            return

        if component_id in self.previous_component_ids:
            return

        self.previous_component_ids.append(component_id)

    @property
    def unmapped_properties(self) -> List[PropertyValidationError]:
        result = []
        for prop in self.input_schema.unmapped_properties:
            result.append(PropertyValidationError(
                msg='Unmapped property',
                loc=[prop.title]  # type: ignore
            ))
        return result

    def set_input_from_previous(self, components: List['Component']) -> bool:
        """
        returns: True if input was changed
        """

        for component in components:
            if component.id not in self.previous_component_ids:
                raise ValueError('Component id not found in previous component ids')

        changed = False

        for prop in self.input_schema.mapped_properties:
            if prop.default:
                path = prop.target_path

                if not prop.target_path:
                    continue

                current_level = self.input_value_dict
                for key in path[:-1]:
                    if key not in current_level:
                        current_level[key] = {}
                    current_level = current_level[key]

                if current_level.get(path[-1]) != prop.default:
                    changed = True

                current_level[path[-1]] = prop.default

        for prev_component in components:
            for prop in self.schema.mapped_properties:
                if prop.source_component_id == prev_component.id:
                    current_level = prev_component.output_value_dict

                    # Find output parameter from output of previous component

                    path = prop.path_from
                    for key in path[:-1]:
                        if key not in current_level:
                            current_level[key] = {}
                        current_level = current_level[key]
                    input_parameter = current_level[path[-1]]

                    # Find and set input parameter for self function

                    path = prop.target_path

                    current_level = self.input_value_dict
                    for key in path[:-1]:
                        if key not in current_level:
                            current_level[key] = {}
                        current_level = current_level[key]

                    if path[-1] not in current_level or current_level[path[-1]] != input_parameter:
                        changed = True

                    current_level[path[-1]] = input_parameter

                    continue

        return changed

    @abstractmethod
    @property
    def jobs_setup_operator_type(self) -> 'JobsSetupOperator':
        ...

    @abstractmethod
    @property
    def job_operator_type(self) -> 'JobOperator':
        ...

    @abstractmethod
    @property
    def output_operator(self) -> 'OutputOperator':
        ...


class JobsSetupOperator(ABC, BaseOperator):
    component_id: uuid.UUID
    repository: WorkflowRepository
    input_changed: bool = False
    '''Whether input was changed after last execution'''

    @apply_defaults
    def __init__(self, component_id: uuid.UUID, task_id: str, **kwargs):
        super().__init__(task_id, **kwargs)

        self.component_id = component_id
        self.repository = WorkflowRepository()

    def pre_execute(self, context: Any):
        component = self.repository.fetch_component(self.component_id)

        prev_components: List[Component] = []

        for previous_component_id in component.previous_component_ids:
            previous_component = self.repository.fetch_component(previous_component_id)

            errors = previous_component.output_errors()
            if errors:
                raise ValueError(errors[0].msg)

            prev_components.append(previous_component)

        self.input_changed = component.set_input_from_previous(prev_components)

        errors = component.output_errors()
        if errors:
            raise ValueError(errors[0].msg)

        self.repository.save_component(component)

    @abstractmethod
    def execute(self, context: Context) -> List[JobId]:
        """
        Setups jobs
        Returns list of job ids
        """
        ...


class JobOperator(ABC, BaseOperator):
    component_id: uuid.UUID
    job_id: JobId

    @apply_defaults
    def __init__(self, job_id: JobId, component_id: uuid.UUID, task_id: str, **kwargs):
        super().__init__(task_id, **kwargs)

        self.component_id = component_id
        self.job_id = job_id

    @abstractmethod
    def execute(self, context: Context) -> Any:
        ...


class OutputOperator(ABC, BaseOperator):
    component_id: uuid.UUID

    @apply_defaults
    def __init__(self, component_id: uuid.UUID, task_id: str, **kwargs):
        super().__init__(task_id, **kwargs)

        self.component_id = component_id

    @abstractmethod
    def execute(self, context: Context) -> Any:
        """
        Post jobs processing
        Setup component output data
        """
        ...
