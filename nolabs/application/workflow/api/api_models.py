import uuid
from enum import Enum
from typing import List, Any, Dict
from uuid import UUID

from pydantic.dataclasses import dataclass


class JobStateEnum(str, Enum):
    RUNNING = 'RUNNING'
    FAILED = 'FAILED'
    COMPLETED = 'COMPLETED'


@dataclass
class GetJobResponse:
    id: uuid.UUID
    state: JobStateEnum
    state_message: str


@dataclass
class GetJobsRequest:
    component_id: UUID


@dataclass
class GetJobsResponse:
    items: List[GetJobResponse]


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
    RUNNING = 'RUNNING'
    FAILED = 'FAILED'
    COMPLETED = 'COMPLETED'


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
    state: ComponentStateEnum
    state_message: str


@dataclass
class GetComponentRequest:
    id: UUID


@dataclass
class AllWorkflowSchemasResponse:
    ids: List[UUID]
