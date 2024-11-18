from __future__ import annotations

__all__ = ["Component"]

import uuid
from abc import abstractmethod
from dataclasses import field
from typing import (
    TYPE_CHECKING,
    Any,
    ClassVar,
    Dict,
    Generic,
    Iterable,
    List,
    Mapping,
    Optional,
    Tuple,
    Type,
    TypeVar,
    Union,
    get_args,
    get_origin,
)
from uuid import UUID

from pydantic import BaseModel, Field, TypeAdapter, ValidationError

from nolabs.domain.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.common import (
    ComponentData,
    PropertyErrorData,
    PropertyValidationError,
)

if TYPE_CHECKING:
    from nolabs.workflow.core.flow import ComponentFlowHandler


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
        return all(
            isinstance(k, key_type) and isinstance(v, value_type)
            for k, v in value.items()
        )
    elif isinstance(value, (list, tuple, set, frozenset)):
        item_type = generic_args[0]

        for item in value:
            if not isinstance(item, item_type):
                try:
                    TypeAdapter(item_type).validate_python(item)
                except ValidationError:
                    return False

        return True
    else:
        return False


class Parameter(BaseModel):
    defs: Optional[Dict[str, "Parameter"]] = Field(alias="$defs", default_factory=dict)
    description: Optional[str] = None
    properties: Dict[str, "Property"] = Field(default_factory=dict)
    required: List[str] = field(default_factory=list)
    title: Optional[str] = None
    type: Optional[Union[str, List[str]]] = None
    anyOf: List[Union["Property", dict]] = Field(default_factory=list)
    default: Optional[Any] = None
    items: Optional[Union["Items", List["Items"]]] = None
    additionalProperties: Optional[Union[bool, "Parameter"]] = True
    format: Optional[str] = None
    const: Optional[Any] = None
    example: Optional[Any] = None

    @classmethod
    def _find_property(
        cls, schema: "Parameter", target_path: List[str]
    ) -> Optional["Property"]:
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
                        return cls._find_property(
                            schema=ref_schema, target_path=target_path
                        )
                    # Property anyOf is not None and can find property type schema
                    if property.anyOf:
                        for any_of_type in property.anyOf:
                            ref = (
                                any_of_type.ref
                                if isinstance(any_of_type, Property)
                                else any_of_type["$ref"]
                            )
                            ref_type_name = cls.get_ref_type_name(ref)
                            if not schema.defs:
                                return None
                            ref_schema = schema.defs[ref_type_name]
                            return cls._find_property(
                                schema=ref_schema, target_path=target_path
                            )
                else:
                    return property

        return None

    @property
    def mapped_properties(self) -> List["Property"]:
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
        return ref.split("/")[-1]

    def try_set_mapping(
        self,
        source_schema: "Parameter",
        component_id: UUID,
        path_from: List[str],
        target_path: List[str],
    ) -> Optional[PropertyValidationError]:
        if not path_from or not target_path:
            raise ValueError("Path from or path to are empty")

        source_property = self._find_property(
            schema=source_schema, target_path=path_from
        )
        if not source_property:
            return PropertyValidationError(
                msg=f"Property does not exist in source schema", loc=path_from
            )

        target_property = self._find_property(schema=self, target_path=target_path)
        if not target_property:
            return PropertyValidationError(
                msg=f"Property does not exist in target schema", loc=target_path
            )

        validation_passed = False

        for any_of in source_property.anyOf:
            if any_of and (
                any_of.type == target_property.type
                or any_of.format == target_property.format
            ):
                validation_passed = True

        if (
            source_property.type == target_property.type
            or source_property.format == target_property.format
        ):
            validation_passed = True

        if not validation_passed:
            return PropertyValidationError(
                msg=f'Properties "{path_from[-1]}" and "{target_path[-1]}" has incompatible types or formats',
                loc=target_path,
            )

        target_property.map(
            source_component_id=component_id,
            path_from=path_from,
            target_path=target_path,
        )

        return None

    def try_set_default(
        self,
        input_type: Type[BaseModel],
        target_path: List[str],
        value: Optional[Any] = None,
    ) -> Optional[PropertyValidationError]:
        if not target_path:
            raise ValueError("Path from or path two are empty")

        target_property = self._find_property(schema=self, target_path=target_path)
        if not target_property:
            return PropertyValidationError(
                msg=f"Property does not exist in target schema", loc=target_path
            )

        annotations = input_type.__annotations__
        current_type = Type

        for path in target_path:
            if path not in annotations:
                raise ValueError("Path is missing from the annotations")

            current_type = annotations[path]
            if hasattr(current_type, "__annotations__"):
                annotations = current_type.__annotations__

        if not is_assignable_to_generic(value, current_type):
            return PropertyValidationError(
                msg=f"Property has incompatible type", loc=target_path
            )

        target_property.default = value
        target_property.target_path = target_path

    @staticmethod
    def get_instance(cls: Type) -> "Parameter":
        if not issubclass(cls, BaseModel):
            raise ValueError(f"Schema must be a subclass of {BaseModel}")

        schema = cls.model_json_schema()

        return Parameter(**schema)

    @property
    def unmapped_properties(self) -> List["Property"]:
        result: List[Property] = []

        if not self.properties:
            return result

        for name, property in self.properties.items():
            if name in self.required and not (
                property.default or property.source_component_id
            ):
                result.append(property)

        if not self.defs:
            return result

        for _, schema in self.defs.items():
            if not schema or not schema.properties:
                continue
            for name, property in schema.properties.items():
                # TODO check default value
                if name in self.required and not (
                    property.default or property.source_component_id
                ):
                    result.append(property)

        return result

    def _validate_dictionary(
        self, t: Type, dictionary: Dict[str, Any]
    ) -> List[PropertyValidationError]:
        try:
            _ = t(**dictionary)
        except ValidationError as e:
            return [
                PropertyValidationError(
                    msg=error["msg"], loc=error["loc"]  # type: ignore
                )
                for error in e.errors()
            ]

        return []


class Property(BaseModel):
    type: Optional[Union[str, List[str]]] = None
    properties: Optional[Dict[str, "Property"]] = Field(default_factory=dict)
    items: Optional[Union["Items", List["Items"]]] = None
    required: List[str] = Field(default_factory=list)
    description: Optional[str] = None
    enum: List[Any] = Field(default_factory=list)
    const: Optional[Any] = None
    format: Optional[str] = None
    default: Optional[Any] = None
    example: Optional[Any] = None
    title: Optional[str] = None
    target_path: List[str] = Field(default_factory=list)
    anyOf: List[Union["Property", dict]] = Field(default_factory=list)
    ref: Optional[str] = Field(alias="$ref", default=None)

    source_component_id: Optional[UUID] = None
    path_from: List[str] = Field(default_factory=list)

    def map(
        self, source_component_id: UUID, path_from: List[str], target_path: List[str]
    ):
        self.source_component_id = source_component_id
        self.path_from = path_from
        self.target_path = target_path


class Items(BaseModel):
    type: Optional[Union[str, List[str]]] = None
    properties: Optional[Dict[str, Property]] = Field(default_factory=dict)
    required: List[str] = Field(default_factory=list)
    items: Optional[Union["Items", List["Items"]]] = None
    description: Optional[str] = None
    enum: List[Any] = Field(default_factory=list)
    ref: Any = Field(alias="$ref", default=None)
    const: Optional[Any] = None
    format: Optional[str] = None
    default: Optional[Any] = None
    example: Optional[Any] = None


TInput = TypeVar("TInput", bound=BaseModel)
TOutput = TypeVar("TOutput", bound=BaseModel)


def is_pydantic_type(t: Any) -> bool:
    return issubclass(type(t), BaseModel) or "__pydantic_post_init__" in t.__dict__


class Component(Generic[TInput, TOutput]):
    id: uuid.UUID

    # region schemas

    output_schema: Parameter
    output_value_dict: Dict[str, Any]

    input_schema: Parameter
    input_value_dict: Dict[str, Any]

    previous_component_ids: List[uuid.UUID] = []

    # endregion

    name: ClassVar[str]
    description: ClassVar[str]

    meta = {"collection": "components"}

    def __init__(
        self,
        id: uuid.UUID,
        input_schema: Optional[Union[Parameter, Dict[str, Any]]] = None,
        output_schema: Optional[Union[Parameter, Dict[str, Any]]] = None,
        input_value_dict: Optional[Dict[str, Any]] = None,
        output_value_dict: Optional[Dict[str, Any]] = None,
        previous_component_ids: Optional[List[uuid]] = None,
    ):
        self.id = id

        if isinstance(input_schema, Mapping):
            self.input_schema = Parameter(**input_schema)
        else:
            self.input_schema = input_schema or Parameter.get_instance(
                cls=self.input_parameter_type
            )

        if isinstance(output_schema, Mapping):
            self.output_schema = Parameter(**output_schema)
        else:
            self.output_schema = output_schema or Parameter.get_instance(
                cls=self.output_parameter_type
            )

        self.input_value_dict = input_value_dict or {}
        self.output_value_dict = output_value_dict or {}
        self.previous_component_ids = previous_component_ids or []

        if not self.input_value_dict and not self.input_errors() and self.input_value:
            self.input_value_dict = self.input_value.model_dump()

        if not self.output_value_dict and not self.output_errors() and self.output_value:
            self.output_value_dict = self.output_value.model_dump()

    @property
    def output_value(self) -> TOutput:
        return self.output_parameter_type(**self.output_value_dict)

    @output_value.setter
    def output_value(self, value: Union[TOutput, Dict[str, Any]]):
        if isinstance(value, dict):
            self.output_value_dict = value
        else:
            self.output_value_dict = value.model_dump()

    def output_errors(self, rse=False) -> List[PropertyValidationError]:
        errors = self.output_schema._validate_dictionary(
            t=self.output_parameter_type, dictionary=self.output_value_dict
        )
        if rse and errors:
            messages = []
            for error in errors:
                message = f'{error.msg}: '
                message += ', '.join(x for x in error.loc)
                messages.append(message)
            raise NoLabsException(ErrorCodes.component_input_invalid, message='; '.join(messages))
        return errors

    @property
    @abstractmethod
    def input_parameter_type(self) -> Type[TInput]: ...

    @property
    @abstractmethod
    def output_parameter_type(self) -> Type[TOutput]: ...

    @property
    def input_value(self) -> TInput:
        return self.input_parameter_type(**self.input_value_dict)

    def input_errors(self, rse=False) -> List[PropertyValidationError]:
        errors = self.input_schema._validate_dictionary(
            t=self.input_parameter_type, dictionary=self.input_value_dict
        )
        if rse and errors:
            messages = []
            for error in errors:
                message = f'Msg: {error.msg}'
                message += ', '.join(x for x in error.loc)
                messages.append(message)
            raise NoLabsException(ErrorCodes.component_input_invalid, message='; '.join(messages))
        return errors

    def try_map_property(
        self, component: "Component", path_from: List[str], target_path: List[str]
    ) -> Optional[PropertyValidationError]:
        return self.input_schema.try_set_mapping(
            source_schema=component.output_schema,
            component_id=component.id,
            path_from=path_from,
            target_path=target_path,
        )

    def try_set_default(
        self, target_path: List[str], value: Any
    ) -> Optional[PropertyValidationError]:
        error = self.input_schema.try_set_default(
            target_path=target_path, value=value, input_type=self.input_parameter_type
        )

        if not error:
            path = target_path

            current_level = self.input_value_dict
            for key in path[:-1]:
                if key not in current_level:
                    current_level[key] = {}
                current_level = current_level[key]

            current_level[path[-1]] = value

        return error

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
            result.append(
                PropertyValidationError(
                    msg="Unmapped property", loc=[prop.title]  # type: ignore
                )
            )
        return result

    def set_input_from_previous(self, components: List["Component"]) -> bool:
        """
        returns: True if input was changed
        """

        for component in components:
            if component.id not in self.previous_component_ids:
                raise ValueError("Component id not found in previous component ids")

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

                    if (
                        path[-1] not in current_level
                        or current_level[path[-1]] != input_parameter
                    ):
                        changed = True

                    current_level[path[-1]] = input_parameter

                    continue

        input_errors = self.input_errors()
        if input_errors:
            message = str([(e.msg, e.loc) for e in input_errors])
            raise NoLabsException(ErrorCodes.component_input_invalid, message=message)

        return changed

    @property
    @abstractmethod
    def component_flow_type(self) -> Type["ComponentFlowHandler"]: ...

    @classmethod
    def restore(cls, data: ComponentData) -> "Component":
        component = ComponentTypeFactory.get_type(data.name)(
            id=data.id,
            input_schema=Parameter(**data.input_schema),
            output_schema=Parameter(**data.output_schema),
            input_value_dict=data.input_value_dict,
            output_value_dict=data.output_value_dict,
            previous_component_ids=data.previous_component_ids,
        )

        return component

    def dump(self, data: ComponentData):
        data.name = self.name
        data.input_schema = self.input_schema.model_dump()
        data.output_schema = self.output_schema.model_dump()
        data.input_value_dict = self.input_value_dict
        data.output_value_dict = self.output_value_dict
        data.previous_component_ids = self.previous_component_ids
        data.input_errors = [
            PropertyErrorData.create(loc=err.loc, msg=err.msg)
            for err in self.input_errors()
        ]
        data.output_errors = [
            PropertyErrorData.create(loc=err.loc, msg=err.msg)
            for err in self.output_errors()
        ]


class ComponentTypeFactory:
    _types: ClassVar[Dict[str, Type[Component]]] = {}

    @classmethod
    def enumerate(cls) -> Iterable[Tuple[str, Type[Component]]]:
        return cls._types.items()

    @classmethod
    def add_type(cls, t: Type[Component]):
        cls._types[t.name] = t

    @classmethod
    def get_type(cls, name: str) -> Type[Component]:
        if not cls._types:
            raise ValueError(
                "You must initialize type factory before working with application"
            )

        if name not in cls._types:
            raise ValueError(f"Cannot find component with name {name}")

        return cls._types[name]

    @classmethod
    def clear(cls):
        cls._types = {}
