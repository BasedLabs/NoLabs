from uuid import UUID

from pydantic.dataclasses import dataclass


@dataclass
class UpdateJobRequest:
    job_name: str


@dataclass
class GetJobMetadataResponse:
    job_id: UUID
    job_name: str
    type: str
