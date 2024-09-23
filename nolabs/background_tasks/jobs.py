import uuid

from prefect import flow, get_client, task
from prefect.exceptions import ObjectNotFound

from nolabs.domain.models.common import Job
from domain.models import ComponentData
from nolabs.infrastructure.log import logger
from nolabs.infrastructure.mongo_connector import mongo_connect
from nolabs.infrastructure.settings import settings


@task(task_run_name="cleanup-orphan-task-run-id")
async def _cleanup_orhpan_task_run_id(task_run_id: uuid.UUID) -> bool:
    try:
        async with get_client() as client:
            _ = await client.read_task_run(task_run_id=task_run_id)
        return True
    except ObjectNotFound:
        return False


@flow(flow_run_name="cleanup-orphan-task-run-ids")
async def cleanup_orhpan_task_run_ids():
    mongo_connect(settings.connection_string)

    for job in Job.objects(task_run_id__ne=None).only("id", "task_run_id"):
        exists = await _cleanup_orhpan_task_run_id(job.task_run_id)
        if not exists:
            logger.info(
                "Task run not found for job, cleaning up",
                extra={"job_id": job.id, "task_run_id": job.task_run_id},
            )
            Job.objects(id=job.id).update(set__task_run_id=None)


@task(task_run_name="cleanup-orphan-flow-run-id")
async def _cleanup_orhpan_flow_run_id(flow_run_id: uuid.UUID) -> bool:
    try:
        async with get_client() as client:
            _ = await client.read_flow_run(flow_run_id=flow_run_id)
        return True
    except ObjectNotFound:
        return False


@flow(flow_run_name="cleanup-orphan-flow-run-ids")
async def cleanup_orhpan_flow_run_ids():
    mongo_connect(settings.connection_string)

    cd: ComponentData
    for cd in ComponentData.objects(flow_run_id__ne=None).only("id", "flow_run_id"):
        exists = await _cleanup_orhpan_task_run_id(cd.flow_run_id)
        if not exists:
            logger.info(
                "Flow run not found for component, cleaning up",
                extra={"component_id": cd.id, "flow_run_id": cd.flow_run_id},
            )
            ComponentData.objects(id=cd.id).update(set__flow_run_id=None)
