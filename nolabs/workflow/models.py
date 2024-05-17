import uuid
from typing import Optional, List, Any, Dict

from pydantic import BaseModel, Field


class Parameter(BaseModel):
    type: str
    format: str
    default: Optional[Any] = None


class Component(BaseModel):
    id: str
    type: str
    input: Dict[str, Parameter]
    output: Dict[str, Parameter]


class Mapping(BaseModel):
    source_path: List[str]
    target_path: List[str]
    source_component_id: uuid.UUID
    error: Optional[str] = None


class WorkflowComponent(BaseModel):
    function_id: uuid.UUID
    title: str
    component_id: str
    job_ids: List[uuid.UUID] = Field(default_factory=list)
    mappings: List[Mapping] = Field(default_factory=list)
    error: Optional[str] = None


class WorkflowSchema(BaseModel):
    experiment_id: uuid.UUID
    error: Optional[str]
    components: List[Component]
    workflow: List[WorkflowComponent]
