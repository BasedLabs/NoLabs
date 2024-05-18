from dataclasses import field
from typing import Optional, Union, Dict, List, Any, Tuple

from pydantic import Field, BaseModel
from pydantic.dataclasses import dataclass

from nolabs.workflow.exceptions import WorkflowException


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

    def __post_model_init__(self, __context):
        def _(property: Property, prop_name):
            property.path_to.append(prop_name)

            if not property.properties:
                return

            for nested_name, nested in property.properties.items():
                _(property=nested, prop_name=nested_name)

        for name, property in self.properties.items():
            _(property=property, prop_name=name)

    @classmethod
    def _find_property(cls, schema: 'ParameterSchema', path_to: List[str]) -> Optional['Property']:
        if not path_to:
            return None

        if not schema.properties:
            return None

        path = path_to[0]
        path_to = path_to[1:]

        for name, property in schema.properties.items():
            # We found property in path
            if path == name:
                # Path to is not empty and we must go deeper
                if path_to:
                    if property.ref:
                        ref_type_name = cls._get_ref_type_name(property.ref)
                        if not schema.defs:
                            return None
                        ref_schema = schema.defs[ref_type_name]
                        return cls._find_property(schema=ref_schema, path_to=path_to)
                    # Property anyOf is not None and can find property type definition
                    if property.anyOf:
                        for any_of_type in property.anyOf:
                            ref = any_of_type.ref if isinstance(any_of_type, Property) else any_of_type['$ref']
                            ref_type_name = cls._get_ref_type_name(ref)
                            if not schema.defs:
                                return None
                            ref_schema = schema.defs[ref_type_name]
                            return cls._find_property(schema=ref_schema, path_to=path_to)
                else:
                    return property

        return None

    @property
    def mapped_properties(self) -> List['Property']:
        result: List[Property] = []

        if not self.properties:
            return result

        for _, property in self.properties.items():
            if property.source_function_id:
                result.append(property)

        if not self.defs:
            return result

        for _, definition in self.defs.items():
            if not definition or not definition.properties:
                continue

            for _, property in definition.properties.items():
                if property.source_function_id:
                    result.append(property)

        return result

    @staticmethod
    def _get_ref_type_name(ref: str) -> str:
        return ref.split('/')[-1]

    def find_property(self, path: List[str]) -> Optional['Property']:
        return self._find_property(schema=self, path_to=path)

    def try_set_mapping(self, source_schema: 'ParameterSchema',
                        function_id: str,
                        path_from: List[str],
                        path_to: List[str]) -> Optional[PropertyValidationError]:
        if not path_from or not path_to:
            raise WorkflowException(
                msg='Path from or path two are empty'
            )

        source_property = self._find_property(schema=source_schema, path_to=path_from)
        if not source_property:
            return PropertyValidationError(
                msg=f'Property does not exist in source schema',
                loc=path_from
            )

        target_property = self._find_property(schema=self, path_to=path_to)
        if not target_property:
            return PropertyValidationError(
                msg=f'Property does not exist in target schema',
                loc=path_to
            )

        validation_passed = False

        for any_of in source_property.anyOf:
            if any_of and (any_of.type == target_property.type or any_of.format == target_property.format):
                validation_passed = True

        if source_property.type == target_property.type or source_property.format == target_property.format:
            validation_passed = True

        if not validation_passed:
            return PropertyValidationError(
                msg=f'Properties "{path_from[-1]}" and "{path_to[-1]}" has incompatible types or formats',
                loc=path_to
            )

        target_property.map(
            source_function_id=function_id,
            path_from=path_from
        )

        return None

    @staticmethod
    def get_instance(cls) -> 'ParameterSchema':
        if not issubclass(cls, BaseModel):
            raise WorkflowException(
                f'Schema must be a subclass of {BaseModel}'
            )

        return ParameterSchema(**cls.schema())

    @property
    def unmapped_properties(self) -> List['Property']:
        result: List[Property] = []

        if not self.properties:
            return result

        for name, property in self.properties.items():
            if name in self.required and not property.source_function_id:
                result.append(property)

        if not self.defs:
            return result

        for _, definition in self.defs.items():
            if not definition or not definition.properties:
                continue
            for name, property in definition.properties.items():
                # TODO check default value
                if name in self.required and not property.source_function_id:
                    result.append(property)

        return result

    def unmap(self, property: Optional['Property'] = None, path_to: Optional[List[str]] = None):
        if not property and not path_to:
            raise WorkflowException(
                msg='You must provide either property or path_to'
            )

        if path_to:
            property = self._find_property(schema=self, path_to=path_to)
            if property:
                property.unmap()
                return

        if property:
            property.unmap()


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
    path_to: List[str] = Field(default_factory=list)
    anyOf: List[Union['Property', dict]] = Field(default_factory=list)
    ref: Optional[str] = Field(alias='$ref', default=None)

    source_function_id: Optional[str] = None
    path_from: List[str] = Field(default_factory=list)

    def unmap(self):
        self.source_function_id = None
        self.path_from = []

    def map(self, source_function_id: str, path_from: List[str]):
        self.source_function_id = source_function_id
        self.path_from = path_from


class Items(BaseModel):
    type: Optional[Union[str, List[str]]] = None
    properties: Optional[Dict[str, Property]] = Field(default_factory=dict)
    required: List[str] = Field(default_factory=list)
    items: Optional[Union['Items', List['Items']]] = None
    description: Optional[str] = None
    enum: List[Any] = Field(default_factory=list)
    const: Optional[Any] = None
    format: Optional[str] = None
    default: Optional[Any] = None
    example: Optional[Any] = None
