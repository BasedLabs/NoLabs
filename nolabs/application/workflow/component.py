__all__ = [
    'Component'
]

import uuid
from abc import abstractmethod, ABC
from dataclasses import field
from typing import Iterable, Tuple
from typing import Optional, Union, Dict, List, Any, Type, get_origin, get_args
from typing import TypeVar, Generic, ClassVar, Mapping
from uuid import UUID

from airflow.models import BaseOperator
from pydantic import Field, BaseModel, ValidationError, parse_obj_as
from pydantic.dataclasses import dataclass

from nolabs.application.workflow.states import ComponentState


def is_assignable_to_generic(value, generic_type):
    origin_type = get_origin(generic_type)

    if not origin_type:
        return isinstance(value, generic_type)

    if not isinstance(value, origin_type):
        return False

    generic_args = get_args(generic_type)
    if not generic_args:
        return True  # No arguments to check, so it's a match

    if isinstance(value, dict):
        key_type, value_type = generic_args
        return all(isinstance(k, key_type) and isinstance(v, value_type) for k, v in value.items())
    elif isinstance(value, (list, tuple, set, frozenset)):
        item_type = generic_args[0]

        for item in value:
            if not isinstance(item, item_type):
                try:
                    parse_obj_as(item_type, item)
                except ValidationError:
                    return False

        return True
    else:
        return False


@dataclass
class PropertyValidationError:
    msg: str
    loc: List[str]

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'{self.msg}: {self.loc}'


class ParameterSchema(BaseModel):
    defs: Optional[Dict[str, 'ParameterSchema']] = Field(alias='$defs', default_factory=dict)
    description: Optional[str] = None
    properties: Dict[str, 'Property'] = Field(default_factory=dict)
    required: List[str] = field(default_factory=list)
    title: Optional[str] = None
    type: Optional[Union[str, List[str]]] = None
    anyOf: List[Union['Property', dict]] = Field(default_factory=list)
    default: Optional[Any] = None
    items: Optional[Union['Items', List['Items']]] = None
    additionalProperties: Optional[Union[bool, 'ParameterSchema']] = True
    format: Optional[str] = None
    const: Optional[Any] = None
    example: Optional[Any] = None

    @classmethod
    def _find_property(cls, schema: 'ParameterSchema', target_path: List[str]) -> Optional['Property']:
        if not target_path:
            return None

        if not schema.properties:
            return None

        path = target_path[0]
        target_path = target_path[1:]

        for name, property in schema.properties.items():
            # We found property in path
            if path == name:
                # Path to is not empty and we must go deeper
                if target_path:
                    if property.ref:
                        ref_type_name = cls.get_ref_type_name(property.ref)
                        if not schema.defs:
                            return None
                        ref_schema = schema.defs[ref_type_name]
                        return cls._find_property(schema=ref_schema, target_path=target_path)
                    # Property anyOf is not None and can find property type schema
                    if property.anyOf:
                        for any_of_type in property.anyOf:
                            ref = any_of_type.ref if isinstance(any_of_type, Property) else any_of_type['$ref']
                            ref_type_name = cls.get_ref_type_name(ref)
                            if not schema.defs:
                                return None
                            ref_schema = schema.defs[ref_type_name]
                            return cls._find_property(schema=ref_schema, target_path=target_path)
                else:
                    return property

        return None

    @property
    def mapped_properties(self) -> List['Property']:
        result: List[Property] = []

        if not self.properties:
            return result

        for _, property in self.properties.items():
            if property.source_component_id or property.default:
                result.append(property)

        if not self.defs:
            return result

        for _, schema in self.defs.items():
            if not schema or not schema.properties:
                continue

            for _, property in schema.properties.items():
                if property.source_component_id or property.default:
                    result.append(property)

        return result

    @staticmethod
    def get_ref_type_name(ref: str) -> str:
        return ref.split('/')[-1]

    def find_property(self, path: List[str]) -> Optional['Property']:
        return self._find_property(schema=self, target_path=path)

    def try_set_mapping(self, source_schema: 'ParameterSchema',
                        component_id: UUID,
                        path_from: List[str],
                        target_path: List[str]) -> Optional[PropertyValidationError]:
        if not path_from or not target_path:
            raise ValueError('Path from or path to are empty')

        source_property = self._find_property(schema=source_schema, target_path=path_from)
        if not source_property:
            return PropertyValidationError(
                msg=f'Property does not exist in source schema',
                loc=path_from
            )

        target_property = self._find_property(schema=self, target_path=target_path)
        if not target_property:
            return PropertyValidationError(
                msg=f'Property does not exist in target schema',
                loc=target_path
            )

        validation_passed = False

        for any_of in source_property.anyOf:
            if any_of and (any_of.type == target_property.type or any_of.format == target_property.format):
                validation_passed = True

        if source_property.type == target_property.type or source_property.format == target_property.format:
            validation_passed = True

        if not validation_passed:
            return PropertyValidationError(
                msg=f'Properties "{path_from[-1]}" and "{target_path[-1]}" has incompatible types or formats',
                loc=target_path
            )

        target_property.map(
            source_component_id=component_id,
            path_from=path_from,
            target_path=target_path
        )

        return None

    def try_set_default(self,
                        input_type: Type[BaseModel],
                        target_path: List[str],
                        value: Optional[Any] = None) -> Optional[PropertyValidationError]:
        if not target_path:
            raise ValueError('Path from or path two are empty')

        target_property = self._find_property(schema=self, target_path=target_path)
        if not target_property:
            return PropertyValidationError(
                msg=f'Property does not exist in target schema',
                loc=target_path
            )

        annotations = input_type.__annotations__
        current_type = Type

        for path in target_path:
            if path not in annotations:
                raise ValueError('Path is missing from the annotations')

            current_type = annotations[path]
            if hasattr(current_type, '__annotations__'):
                annotations = current_type.__annotations__

        if not is_assignable_to_generic(value, current_type):
            return PropertyValidationError(
                msg=f'Property has incompatible type',
                loc=target_path
            )

        target_property.default = value
        target_property.target_path = target_path

    @staticmethod
    def get_instance(cls: Type) -> 'ParameterSchema':
        if not issubclass(cls, BaseModel):
            raise ValueError(
                f'Schema must be a subclass of {BaseModel}'
            )

        schema = cls.schema()

        return ParameterSchema(**schema)

    @property
    def unmapped_properties(self) -> List['Property']:
        result: List[Property] = []

        if not self.properties:
            return result

        for name, property in self.properties.items():
            if name in self.required and not (property.default or property.source_component_id):
                result.append(property)

        if not self.defs:
            return result

        for _, schema in self.defs.items():
            if not schema or not schema.properties:
                continue
            for name, property in schema.properties.items():
                # TODO check default value
                if name in self.required and not (property.default or property.source_component_id):
                    result.append(property)

        return result

    def unmap(self, property: Optional['Property'] = None, target_path: Optional[List[str]] = None):
        if not property and not target_path:
            raise ValueError(
                'You must provide either property or target_path'
            )

        if target_path:
            property = self._find_property(schema=self, target_path=target_path)
            if property:
                property.unmap()
                return

        if property:
            property.unmap()

    def validate_dictionary(self, t: Type, dictionary: Dict[str, Any]) -> List[PropertyValidationError]:
        try:
            _ = t(**dictionary)
        except ValidationError as e:
            return [
                PropertyValidationError(
                    msg=error['msg'],
                    loc=error['loc']  # type: ignore
                )
                for error in e.errors()
            ]

        return []


class Property(BaseModel):
    type: Optional[Union[str, List[str]]] = None
    properties: Optional[Dict[str, 'Property']] = Field(default_factory=dict)
    items: Optional[Union['Items', List['Items']]] = None
    required: List[str] = Field(default_factory=list)
    description: Optional[str] = None
    enum: List[Any] = Field(default_factory=list)
    const: Optional[Any] = None
    format: Optional[str] = None
    default: Optional[Any] = None
    example: Optional[Any] = None
    title: Optional[str] = None
    target_path: List[str] = Field(default_factory=list)
    anyOf: List[Union['Property', dict]] = Field(default_factory=list)
    ref: Optional[str] = Field(alias='$ref', default=None)

    source_component_id: Optional[UUID] = None
    path_from: List[str] = Field(default_factory=list)

    def unmap(self):
        self.source_component_id = None
        self.path_from = []

    def map(self, source_component_id: UUID, path_from: List[str], target_path: List[str]):
        self.source_component_id = source_component_id
        self.path_from = path_from
        self.target_path = target_path


class Items(BaseModel):
    type: Optional[Union[str, List[str]]] = None
    properties: Optional[Dict[str, Property]] = Field(default_factory=dict)
    required: List[str] = Field(default_factory=list)
    items: Optional[Union['Items', List['Items']]] = None
    description: Optional[str] = None
    enum: List[Any] = Field(default_factory=list)
    ref: Any = Field(alias='$ref', default=None)
    const: Optional[Any] = None
    format: Optional[str] = None
    default: Optional[Any] = None
    example: Optional[Any] = None


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
            self.input_schema = input_schema or ParameterSchema.get_instance(cls=self.input_parameter_type)

        if isinstance(output_schema, Mapping):
            self.output_schema = ParameterSchema(**output_schema)
        else:
            self.output_schema = output_schema or ParameterSchema.get_instance(cls=self.output_parameter_type)

        self.input_value_dict = input_value_dict or {}
        self.output_value_dict = output_value_dict or {}
        self.previous_component_ids = previous_component_ids or []

    @property
    def output_value(self) -> TOutput:
        return self.output_parameter_type(**self.output_value_dict)

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

        component = ComponentTypeFactory.get_type(state.name)(
            id=state.id,
            experiment_id=state.experiment_id,
            job_ids=state.job_ids,
            input_schema=ParameterSchema(**state.input_schema),
            output_schema=ParameterSchema(**state.output_schema),
            input_value_dict=state.input_value_dict,
            output_value_dict=state.output_value_dict,
            previous_component_ids=state.previous_component_ids
        )

        component._state = state
        return component

    def save(self):
        self._state.set_component(component=self)
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
