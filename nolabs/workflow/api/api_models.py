import uuid
from enum import Enum
from typing import Any, Dict, List, Optional
from uuid import UUID

from pydantic.dataclasses import dataclass


class JobStateEnum(str, Enum):
    RUNNING = "RUNNING"
    FAILED = "FAILED"
    COMPLETED = "COMPLETED"


@dataclass
class GetJobState:
    id: uuid.UUID
    component_id: uuid.UUID
    state: Optional[JobStateEnum]
    state_message: Optional[str]


@dataclass
class GetJobRequest:
    job_id: uuid.UUID


@dataclass
class PropertyErrorResponse:
    loc: List[str]
    msg: str


@dataclass
class ResetWorkflowRequest:
    workflow_id: UUID


@dataclass
class StartWorkflowComponentRequest:
    workflow_id: UUID
    component_id: UUID


class ComponentStateEnum(str, Enum):
    RUNNING = "RUNNING"
    FAILED = "FAILED"
    COMPLETED = "COMPLETED"


@dataclass
class GetComponentResponse:
    id: UUID
    input_schema: Dict[str, Any]
    output_schema: Dict[str, Any]
    input_value_dict: Dict[str, Any]
    output_value_dict: Dict[str, Any]
    previous_component_ids: List[UUID]
    input_errors: List[PropertyErrorResponse]
    output_errors: List[PropertyErrorResponse]
    state: Optional[ComponentStateEnum]
    state_message: Optional[str]
    job_ids: List[uuid.UUID]


@dataclass
class GetComponentRequest:
    id: UUID


@dataclass
class AllWorkflowSchemasResponse:
    ids: List[UUID]
