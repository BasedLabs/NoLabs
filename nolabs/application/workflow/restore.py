__all__ = ["restore_statuses"]

from prefect import get_run_logger, task
from prefect.client.orchestration import get_client
from prefect.client.schemas import StateType

from nolabs.workflow.data import ComponentData, JobRunData


@task
async def _restore_statuses():
    logger = get_run_logger()

    async with get_client() as client:
        logger.info("Starting restore task")

        jobs = JobRunData.objects(state=StateType.RUNNING)

        job: JobRunData
        for job in jobs:
            logger.info(f"Restoring job state %s", job.id)
            task_run = await client.read_task_run(job.task_run_id)

            if task_run.state_type.value != job.state:
                job.state = task_run.state_type.value
                job.state_message = task_run.state.message

        components = ComponentData.objects(state=StateType.RUNNING)

        component: ComponentData
        for component in components:
            logger.info(f"Restoring component state %s", job.id)

            flow_run = await client.read_flow_run(component.flow_run_id)

            if flow_run.state_type.value != component.state:
                component.state = flow_run.state_type.value
                component.state_message = flow_run.state.message


def restore_statuses():
    _restore_statuses.delay()
