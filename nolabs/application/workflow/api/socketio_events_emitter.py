from typing import List
from uuid import UUID

from infrastructure import socket_server


def emit_start_job_event(experiment_id: UUID, job_id: UUID):
    socket_server.emit_event(name='job_started', id=str(experiment_id), data={'job_id': job_id})


def emit_finish_job_event(experiment_id: UUID, job_id: UUID):
    socket_server.emit_event(name='job_finished', id=str(experiment_id), data={'job_id': job_id})


def emit_start_component_event(experiment_id: UUID, component_id: UUID):
    socket_server.emit_event(name='component_started', id=str(experiment_id), data={'component_id': component_id})


def emit_finish_component_event(experiment_id: UUID, component_id: UUID):
    socket_server.emit_event(name='component_finished', id=str(experiment_id), data={'component_id': component_id})


def emit_component_jobs_event(experiment_id: UUID, component_id: UUID, job_ids: List[UUID]):
    socket_server.emit_event(name='component_jobs', id=str(experiment_id),
                             data={'component_id': component_id, 'job_ids': job_ids})
