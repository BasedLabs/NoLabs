from dataclasses import field
from typing import Optional, Union, Dict, List, Any
from uuid import UUID

from pydantic import Field, BaseModel
from pydantic.dataclasses import dataclass


@dataclass
class SchemaValidationIssue:
    msg: str
    loc: List[str]


class Schema(BaseModel):
    defs: Optional[Dict[str, 'Schema']] = Field(alias='$defs', default_factory=dict)
    description: Optional[str] = None
    properties: Optional[Dict[str, 'Property']] = Field(default_factory=dict)
    required: List[str] = field(default_factory=list)
    title: Optional[str] = None
    type: Optional[Union[str, List[str]]] = None
    anyOf: List[Union['Property', dict]] = Field(default_factory=list)
    default: Optional[Any] = None
    items: Optional[Union['Items', List['Items']]] = None
    additionalProperties: Optional[Union[bool, 'Schema']] = True
    format: Optional[str] = None
    const: Optional[Any] = None
    example: Optional[Any] = None

    def _find_property(self, schema: 'Schema', path_to: List[str]) -> Optional['Property']:
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
                        ref_type_name = self._get_ref_type_name(property.ref)
                        if not self.defs:
                            return None
                        ref_schema = self.defs[ref_type_name]
                        return self._find_property(schema=ref_schema, path_to=path_to)
                    # Property anyOf is not None and can find property type definition
                    if property.anyOf:
                        for any_of_type in property.anyOf:
                            ref = any_of_type.ref if isinstance(any_of_type, Property) else any_of_type['$ref']
                            ref_type_name = self._get_ref_type_name(ref)
                            if not self.defs:
                                return None
                            ref_schema = self.defs[ref_type_name]
                            return self._find_property(schema=ref_schema, path_to=path_to)
                else:
                    return property

        return None

    @property
    def mapped_properties(self) -> List['Property']:
        result: List[Property] = []

        for _, property in self.properties.items():
            if property.mapping_function_id:
                result.append(property)

        for _, definition in self.defs.items():
            for _, property in definition.properties.items():
                if property.mapping_function_id:
                    result.append(property)

        return result

    def _get_ref_type_name(self, ref: str) -> str:
        return ref.split('/')[-1]

    def find_property(self, path: List[str]) -> Optional['Property']:
        return self._find_property(schema=self, path_to=path)

    def remove_mapping(self, path_to: List[str]):
        for property in self.mapped_properties:
            if property.path_to == path_to:
                property.mapping_function_id = None
                property.path_to = []
                property.path_from = []

    def try_set_mapping(self, source_schema: 'Schema',
                        function_id: UUID,
                        path_from: List[str],
                        path_to: List[str]) -> Optional[SchemaValidationIssue]:
        if not path_from or not path_to:
            raise ValueError('Path from or path two are empty')

        source_property = self._find_property(schema=source_schema, path_to=path_from)
        if not source_property:
            return SchemaValidationIssue(
                msg=f'Property does not exist in source schema',
                loc=path_from
            )
        target_property = self._find_property(schema=self, path_to=path_to)
        if not target_property:
            return SchemaValidationIssue(
                msg=f'Property does not exist in target schema',
                loc=path_to
            )

        if source_property.type != target_property.type or source_property.type != target_property.type:
            return SchemaValidationIssue(
                msg=f'Properties "{path_from[-1]}" and "{path_to[-1]}" has incompatible types or formats',
                loc=path_to
            )

        target_property.mapping_function_id = function_id
        target_property.path_from = path_from
        target_property.path_to = path_to

        return None

    @staticmethod
    def get_schema(cls) -> 'Schema':
        if not issubclass(cls, BaseModel):
            raise ValueError('Schema must be a subclass of BaseModel')

        return Schema(**cls.schema())

    @property
    def required_are_mapped(self) -> bool:



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
    anyOf: List[Union['Property', dict]] = Field(default_factory=list)
    ref: Optional[str] = Field(alias='$ref', default=None)

    mapping_function_id: Optional[UUID] = None
    path_from: List[str] = Field(default_factory=list)
    path_to: List[str] = Field(default_factory=list)


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
