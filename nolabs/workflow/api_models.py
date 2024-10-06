import uuid
from enum import Enum
from typing import Any, Dict, List, Optional
from uuid import UUID

from pydantic import BaseModel
from pydantic.dataclasses import dataclass


class JobStateEnum(str, Enum):
    RUNNING = "RUNNING"
    FAILED = "FAILED"
    COMPLETED = "COMPLETED"
    CANCELLED = "CANCELLED"
    UNKNOWN = "UNKNOWN"


class GetJobState(BaseModel):
    id: uuid.UUID
    component_id: uuid.UUID
    state: Optional[JobStateEnum]
    state_message: Optional[str]


class GetJobRequest(BaseModel):
    job_id: uuid.UUID


class PropertyErrorResponse(BaseModel):
    loc: List[str]
    msg: str


class ResetWorkflowRequest(BaseModel):
    experiment_id: UUID


class StartWorkflowComponentRequest(BaseModel):
    experiment_id: UUID
    component_id: UUID


class ComponentStateEnum(str, Enum):
    RUNNING = "RUNNING"
    FAILED = "FAILED"
    COMPLETED = "COMPLETED"
    CANCELLED = "CANCELLED"
    UNKNOWN = "UNKNOWN"


class GetComponentResponse(BaseModel):
    id: UUID
    input_schema: Dict[str, Any]
    output_schema: Dict[str, Any]
    input_value_dict: Dict[str, Any]
    output_value_dict: Dict[str, Any]
    previous_component_ids: List[UUID]
    input_errors: List[PropertyErrorResponse]
    output_errors: List[PropertyErrorResponse]
    state: ComponentStateEnum
    job_ids: List[uuid.UUID]
    state_message: Optional[str] = None


class GetComponentRequest(BaseModel):
    id: UUID


class AllWorkflowSchemasResponse(BaseModel):
    ids: List[UUID]
