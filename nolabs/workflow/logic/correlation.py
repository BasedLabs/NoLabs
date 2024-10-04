import uuid
from typing import Optional


def make_correlation_id(experiment_id: uuid.UUID, component_id: uuid.UUID, job_id: Optional[uuid.UUID] = None):
    if not job_id:
        return f"{str(experiment_id)}:{str(component_id)}"

    return f"{str(experiment_id)}:{str(component_id)}:{str(job_id)}"


def unpack_correlation_id(correlation_id: str) -> (uuid.UUID, uuid.UUID, Optional[uuid.UUID]):
    parts = correlation_id.split(":")

    if len(parts) == 2:
        return uuid.UUID(parts[0]), uuid.UUID(parts[1]), None

    return uuid.UUID(parts[0]), uuid.UUID(parts[1]), uuid.UUID(parts[2])
