import uuid
from typing import Optional, List, Any, Dict, Union

from pydantic import BaseModel, Field


class PropertyModel(BaseModel):
    type: Optional[Union[str, List[str]]] = None
    properties: Optional[Dict[str, 'PropertyModel']] = Field(default_factory=dict)
    required: List[str] = Field(default_factory=list)
    description: Optional[str] = None
    enum: List[Any] = Field(default_factory=list)
    const: Optional[Any] = None
    format: Optional[str] = None
    default: Optional[Any] = None
    example: Optional[Any] = None
    title: Optional[str] = None
    anyOf: List[Union['PropertyModel', dict]] = Field(default_factory=list)
    ref: Optional[str] = Field(default=None)


class ComponentModel(BaseModel):
    id: str
    name: str
    input: Dict[str, PropertyModel]
    output: Dict[str, PropertyModel]


class MappingModel(BaseModel):
    source_path: List[str]
    target_path: List[str]
    source_component_id: str
    error: Optional[str] = None


class WorkflowComponentModel(BaseModel):
    function_id: uuid.UUID
    title: str
    component_id: str
    job_ids: List[uuid.UUID] = Field(default_factory=list)
    mappings: List[MappingModel] = Field(default_factory=list)
    error: Optional[str] = None


class WorkflowSchemaModel(BaseModel):
    experiment_id: uuid.UUID
    components: List[ComponentModel]
    workflow_components: List[WorkflowComponentModel]
    error: Optional[str] = None
    valid: bool = True
