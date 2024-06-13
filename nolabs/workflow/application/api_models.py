from typing import List, Optional, Any, Dict
from uuid import UUID

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
    jobs_errors: List[JobErrorResponse]
    input_property_errors: List[InputPropertyErrorResponse]
    last_exceptions: List[str]



@dataclass
class GetComponentStateRequest:
    component_id: UUID


@dataclass
class AllWorkflowSchemasResponse:
    ids: List[UUID]