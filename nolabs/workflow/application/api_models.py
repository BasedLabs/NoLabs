from typing import List, Optional, Any, Dict
from uuid import UUID

from pydantic.dataclasses import dataclass


@dataclass
class GetComponentJobIdsRequest:
    component_id: UUID


@dataclass
class GetComponentStateResponse:
    input_dict: Dict[str, Any]
    output_dict: Dict[str, Any]
    job_ids: List[UUID]


@dataclass
class GetComponentStateRequest:
    component_id: UUID


@dataclass
class AllWorkflowSchemasResponse:
    ids: List[UUID]