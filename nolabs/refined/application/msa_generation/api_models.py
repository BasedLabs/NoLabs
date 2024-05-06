__all__ = [
    'PredictMsaRequest',
    'PredictMsaResponse',
    'GetJobStatusResponse'
]

from uuid import UUID

from pydantic.dataclasses import dataclass


@dataclass
class PredictMsaRequest:
    job_id: UUID
    experiment_id: UUID
    protein_id: UUID


@dataclass
class PredictMsaResponse:
    msa: str | None


@dataclass
class GetJobStatusResponse:
    running: bool


