import uuid
from typing import Any, Dict, List, Optional, Union

from pydantic import BaseModel, Field


class ItemsSchema(BaseModel):
    type: Optional[Union[str, List[str]]] = None
    properties: Optional[Dict[str, "PropertySchema"]] = Field(default_factory=dict)
    required: List[str] = Field(default_factory=list)
    description: Optional[str] = None
    enum: List[Any] = Field(default_factory=list)
    ref: Any = Field(alias="$ref", default=None)
    const: Optional[Any] = None
    format: Optional[str] = None
    default: Optional[Any] = None
    example: Optional[Any] = None
    items: Optional[Union["ItemsSchema", List["ItemsSchema"]]] = None


class PropertySchema(BaseModel):
    type: Optional[Union[str, List[str]]] = None
    properties: Optional[Dict[str, "PropertySchema"]] = Field(default_factory=dict)
    required: List[str] = Field(default_factory=list)
    description: Optional[str] = None
    enum: List[Any] = Field(default_factory=list)
    const: Optional[Any] = None
    format: Optional[str] = None
    default: Optional[Any] = None
    example: Optional[Any] = None
    title: Optional[str] = None
    anyOf: List[Union["PropertySchema", dict]] = Field(default_factory=list)
    ref: Optional[str] = Field(default=None)
    items: Optional[Union["ItemsSchema", List["ItemsSchema"]]] = None


class ComponentSchemaTemplate(BaseModel):
    name: str
    input: Dict[str, PropertySchema]
    output: Dict[str, PropertySchema]
    description: Optional[str] = None


class MappingSchema(BaseModel):
    source_path: List[str]
    target_path: List[str]
    source_component_id: uuid.UUID
    error: Optional[str] = None


class DefaultSchema(BaseModel):
    target_path: List[str]
    value: Optional[Any] = None
    error: Optional[str] = None


class ComponentSchema(BaseModel):
    name: str
    component_id: uuid.UUID
    error: Optional[str] = None
    mappings: List[MappingSchema] = Field(default_factory=list)
    defaults: List[DefaultSchema] = Field(default_factory=list)
    x: float = 0.0
    y: float = 0.0


class WorkflowSchema(BaseModel):
    experiment_id: uuid.UUID
    component_templates: List[ComponentSchemaTemplate]
    components: List[ComponentSchema]
    error: Optional[str] = None
    valid: bool = True
