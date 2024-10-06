from typing import List
from uuid import UUID

from nolabs.infrastructure import socket_server
from nolabs.infrastructure.log import logger


def emit_start_job_event(experiment_id: UUID, component_id: UUID, job_id: UUID):
    event_name = "job_started"

    logger.info(
        f"Emitting event",
        extra={
            "event_name": event_name,
            "experiment_id": experiment_id,
            "component_id": component_id,
            "job_id": job_id,
        },
    )
    socket_server.enqueue_event(
        name=event_name,
        room_id=str(experiment_id),
        data={"component_id": str(component_id), "job_id": str(job_id)},
    )


def emit_finish_job_event(experiment_id: UUID, component_id: UUID, job_id: UUID):
    event_name = "job_finished"

    logger.info(
        f"Emitting event",
        extra={
            "event_name": event_name,
            "experiment_id": experiment_id,
            "component_id": component_id,
            "job_id": job_id,
        },
    )

    socket_server.enqueue_event(
        name="job_finished",
        room_id=str(experiment_id),
        data={"component_id": str(component_id), "job_id": str(job_id)},
    )


def emit_start_component_event(experiment_id: UUID, component_id: UUID):
    event_name = "component_started"

    logger.info(
        f"Emitting event",
        extra={
            "event_name": event_name,
            "experiment_id": experiment_id,
            "component_id": component_id,
        },
    )

    socket_server.enqueue_event(
        name="component_started",
        room_id=str(experiment_id),
        data={"component_id": str(component_id)},
    )


def emit_finish_component_event(experiment_id: UUID, component_id: UUID):
    event_name = "component_finished"

    logger.info(
        f"Emitting event",
        extra={
            "event_name": event_name,
            "experiment_id": experiment_id,
            "component_id": component_id,
        },
    )

    socket_server.enqueue_event(
        name="component_finished",
        room_id=str(experiment_id),
        data={"component_id": str(component_id)},
    )


def emit_component_jobs_event(
    experiment_id: UUID, component_id: UUID, job_ids: List[UUID]
):
    event_name = "component_jobs"

    logger.info(
        f"Emitting event",
        extra={
            "event_name": event_name,
            "experiment_id": experiment_id,
            "component_id": component_id,
            "job_ids": job_ids,
        },
    )

    socket_server.enqueue_event(
        name="component_jobs",
        room_id=str(experiment_id),
        data={
            "component_id": str(component_id),
            "job_ids": [str(job_id) for job_id in job_ids],
        },
    )