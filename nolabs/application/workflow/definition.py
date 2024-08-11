import uuid
from typing import Optional, List, Any, Dict, Union

from pydantic import BaseModel, Field


class ItemsTemplate(BaseModel):
    type: Optional[Union[str, List[str]]] = None
    properties: Optional[Dict[str, 'PropertyTemplate']] = Field(default_factory=dict)
    required: List[str] = Field(default_factory=list)
    description: Optional[str] = None
    enum: List[Any] = Field(default_factory=list)
    ref: Any = Field(alias='$ref', default=None)
    const: Optional[Any] = None
    format: Optional[str] = None
    default: Optional[Any] = None
    example: Optional[Any] = None
    items: Optional[Union['ItemsTemplate', List['ItemsTemplate']]] = None


class PropertyTemplate(BaseModel):
    type: Optional[Union[str, List[str]]] = None
    properties: Optional[Dict[str, 'PropertyTemplate']] = Field(default_factory=dict)
    required: List[str] = Field(default_factory=list)
    description: Optional[str] = None
    enum: List[Any] = Field(default_factory=list)
    const: Optional[Any] = None
    format: Optional[str] = None
    default: Optional[Any] = None
    example: Optional[Any] = None
    title: Optional[str] = None
    anyOf: List[Union['PropertyTemplate', dict]] = Field(default_factory=list)
    ref: Optional[str] = Field(default=None)
    items: Optional[Union['ItemsTemplate', List['ItemsTemplate']]] = None


class ComponentTemplate(BaseModel):
    name: str
    input: Dict[str, PropertyTemplate]
    output: Dict[str, PropertyTemplate]
    description: Optional[str] = None


class MappingDefinition(BaseModel):
    source_path: List[str]
    target_path: List[str]
    source_component_id: uuid.UUID
    error: Optional[str] = None


class DefaultDefinition(BaseModel):
    target_path: List[str]
    value: Optional[Any] = None
    error: Optional[str] = None


class ComponentDefinition(BaseModel):
    name: str
    component_id: uuid.UUID
    error: Optional[str] = None
    mappings: List[MappingDefinition] = Field(default_factory=list)
    defaults: List[DefaultDefinition] = Field(default_factory=list)
    x: float = 0.0
    y: float = 0.0


class WorkflowDefinition(BaseModel):
    workflow_id: uuid.UUID
    component_templates: List[ComponentTemplate]
    components: List[ComponentDefinition]
    error: Optional[str] = None
    valid: bool = True
