__all__ = [
    'Component'
]

import uuid
from abc import abstractmethod, ABC
from datetime import datetime
from typing import Dict, Any, List, Iterable, Tuple
from typing import Optional, Type, Union, TypeVar, Generic, ClassVar, Mapping

from airflow.models import BaseOperator
from mongoengine import Document, UUIDField, ReferenceField, CASCADE, DictField, StringField, IntField, \
    ListField, EmbeddedDocumentListField, DateTimeField, EmbeddedDocument
from pydantic import BaseModel
from pydantic.dataclasses import dataclass

from nolabs.application.workflow.api import WorkflowState
from nolabs.application.workflow.properties import ParameterSchema, PropertyValidationError


class InputPropertyErrorDbModel(EmbeddedDocument):
    loc: List[str] = ListField(StringField())
    msg: str = StringField()

    @classmethod
    def create(cls, loc: List[str], msg: str) -> 'InputPropertyErrorDbModel':
        return InputPropertyErrorDbModel(
            loc=loc,
            msg=msg
        )


class ComponentState(Document):
    id: uuid.UUID = UUIDField(primary_key=True)
    experiment_id: uuid.UUID = UUIDField(required=True)
    workflow: WorkflowState = ReferenceField(WorkflowState, reverse_delete_rule=CASCADE)

    input_property_errors: List[InputPropertyErrorDbModel] = EmbeddedDocumentListField(InputPropertyErrorDbModel)

    job_ids: List[uuid.UUID] = ListField(UUIDField())
    last_jobs_count: int = IntField()

    last_executed_at: datetime = DateTimeField()

    name: str = StringField()

    # region component fields

    input_schema: Dict[str, Any] = DictField()
    output_schema: Dict[str, Any] = DictField()
    input_value_dict: Dict[str, Any] = DictField()
    output_value_dict: Dict[str, Any] = DictField()
    previous_component_ids: List[uuid.UUID] = ListField(UUIDField())

    #

    meta = {'collection': 'components'}

    @classmethod
    def create(cls,
               id: uuid.UUID,
               workflow: WorkflowState,
               component: 'Component',
               job_ids: List[uuid.UUID]):
        return ComponentState(
            id=id,
            experiment_id=workflow.experiment.id,
            workflow=workflow,
            input_schema=component.input_schema.dict(),
            output_schema=component.output_schema.dict(),
            input_value_dict=component.input_value_dict,
            output_value_dict=component.output_value_dict,
            previous_component_ids=component.previous_component_ids,
            name=component.name,
            job_ids=job_ids
        )

    def set_component(self, component: 'Component'):
        self.input_schema = component.input_schema.dict()
        self.output_schema = component.output_schema.dict()
        self.input_value_dict = component.input_value_dict
        self.output_value_dict = component.output_value_dict
        self.previous_component_ids = component.previous_component_ids


@dataclass
class JobValidationError:
    job_id: uuid.UUID
    msg: str


TInput = TypeVar('TInput', bound=BaseModel)
TOutput = TypeVar('TOutput', bound=BaseModel)


def is_pydantic_type(t: Any) -> bool:
    return issubclass(type(t), BaseModel) or '__pydantic_post_init__' in t.__dict__


class Component(ABC, Generic[TInput, TOutput]):
    id: uuid.UUID
    experiment_id: uuid.UUID

    job_ids: List[uuid.UUID]

    # region schemas

    output_schema: ParameterSchema
    output_value_dict: Dict[str, Any]

    input_schema: ParameterSchema
    input_value_dict: Dict[str, Any]

    previous_component_ids: List[uuid.UUID] = []

    # endregion

    name: ClassVar[str]
    description: ClassVar[str]

    _state: ComponentState

    def __init__(self,
                 id: uuid.UUID,
                 experiment_id: uuid.UUID,
                 job_ids: Optional[List[uuid.UUID]] = None,
                 input_schema: Optional[Union[ParameterSchema, Dict[str, Any]]] = None,
                 output_schema: Optional[Union[ParameterSchema, Dict[str, Any]]] = None,
                 input_value_dict: Optional[Dict[str, Any]] = None,
                 output_value_dict: Optional[Dict[str, Any]] = None,
                 previous_component_ids: Optional[List[uuid]] = None):
        self.id = id
        self.experiment_id = experiment_id
        self.job_ids = job_ids or []

        if isinstance(input_schema, Mapping):
            self.input_schema = ParameterSchema(**input_schema)
        else:
            self.input_schema = input_schema or ParameterSchema.get_instance(cls=self.input_parameter_type())

        if isinstance(output_schema, Mapping):
            self.output_schema = ParameterSchema(**output_schema)
        else:
            self.output_schema = output_schema or ParameterSchema.get_instance(cls=self.output_parameter_type())

        self.input_value_dict = input_value_dict or {}
        self.output_value_dict = output_value_dict or {}
        self.previous_component_ids = previous_component_ids or []

    @property
    def output_value(self) -> TOutput:
        return self.output_parameter_type()(**self.output_value_dict)

    @output_value.setter
    def output_value(self, value: Union[TOutput, Dict[str, Any]]):
        if value is dict:
            self.output_value_dict = value
        else:
            self.output_value_dict = value.dict()

        if self.output_errors():
            raise ValueError('There were issues while settings the output')

    def output_errors(self) -> List[PropertyValidationError]:
        return self.output_schema.validate_dictionary(t=self.output_parameter_type, dictionary=self.output_value_dict)

    @property
    @abstractmethod
    def input_parameter_type(self) -> Type[TInput]:
        ...

    @property
    @abstractmethod
    def output_parameter_type(self) -> Type[TOutput]:
        ...

    @property
    def input_value(self) -> TInput:
        return self.input_parameter_type()(**self.input_value_dict)

    def input_errors(self):
        return self.input_schema.validate_dictionary(t=self.input_parameter_type(), dictionary=self.input_value_dict)

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
                                                 input_type=self.input_parameter_type())

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
            for prop in self.input_schema.mapped_properties:
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

    @property
    @abstractmethod
    def setup_operator_type(self) -> Type[BaseOperator]:
        ...

    @property
    @abstractmethod
    def job_operator_type(self) -> Optional[Type[BaseOperator]]:
        ...

    @property
    @abstractmethod
    def output_operator_type(self) -> Type[BaseOperator]:
        ...

    @classmethod
    def get(cls, id: uuid.UUID) -> Optional['Component']:
        state: ComponentState = ComponentState.objects.with_id(id)
        if state is None:
            return None
        state: ComponentState = ComponentState.objects.with_id(id)

        if not state:
            return

        return ComponentTypeFactory.get_type(state.name)(
            id=state.id,
            experiment_id=state.experiment_id,
            job_ids=state.job_ids,
            input_schema=ParameterSchema(**state.input_schema),
            output_schema=ParameterSchema(**state.output_schema),
            input_value_dict=state.input_value_dict,
            output_value_dict=state.output_value_dict,
            previous_component_ids=state.previous_component_ids
        )

    def save(self):
        self._state.save()


class ComponentTypeFactory:
    _types: ClassVar[Dict[str, Type[Component]]] = {}

    @classmethod
    def enumerate(cls) -> Iterable[Tuple[str, Type[Component]]]:
        return cls._types.items()

    @classmethod
    def set_types(cls, types: Dict[str, Type[Component]]):
        cls._types = types

    @classmethod
    def get_type(cls, name: str) -> Type[Component]:
        if not cls._types:
            raise ValueError('You must initialize type factory before working with application')

        if name not in cls._types:
            raise ValueError(f"Cannot find component with name {name}")

        return cls._types[name]
