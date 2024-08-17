from typing import List, Any, Dict
from uuid import UUID

from pydantic import Field
from pydantic.dataclasses import dataclass


@dataclass
class GetComponentJobIdsRequest:
    component_id: UUID


@dataclass
class JobErrorResponse:
    msg: str
    job_id: UUID


@dataclass
class InputPropertyErrorResponse:
    loc: List[str]
    msg: str


@dataclass
class ResetWorkflowRequest:
    workflow_id: UUID


@dataclass
class StartWorkflowComponentRequest:
    workflow_id: UUID
    component_id: UUID


@dataclass
class GetComponentStateResponse:
    input_dict: Dict[str, Any]
    output_dict: Dict[str, Any]
    job_ids: List[UUID]
    input_property_errors: List[InputPropertyErrorResponse]
    last_exceptions: List[str]
    jobs_errors: List[JobErrorResponse] = Field(default_factory=list)


@dataclass
class GetComponentStateRequest:
    component_id: UUID


@dataclass
class AllWorkflowDefinitionsResponse:
    ids: List[UUID]
