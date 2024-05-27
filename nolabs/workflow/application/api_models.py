from typing import List, Optional, Any, Dict
from uuid import UUID

from pydantic.dataclasses import dataclass


@dataclass
class GetComponentJobIdsResponse:
    job_ids: List[UUID]


@dataclass
class GetComponentJobIdsRequest:
    name: str
    component_id: UUID


@dataclass
class GetComponentParametersResponse:
    input_dict: Dict[str, Any]
    output_dict: Dict[str, Any]


@dataclass
class GetComponentParametersRequest:
    name: str
    component_id: UUID