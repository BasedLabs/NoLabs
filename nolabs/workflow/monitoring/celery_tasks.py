from celery import Celery
from celery.result import AsyncResult
from celery.states import FAILURE, UNREADY_STATES

from nolabs.infrastructure.celery_app_factory import get_celery_app
from nolabs.infrastructure.redis_client_factory import redlock
from nolabs.infrastructure.settings import settings
from nolabs.workflow.core import Tasks
from nolabs.workflow.monitoring.tasks import OrphanedTasksTracker


def register_monitoring_celery_tasks(celery: Celery):
    @celery.task(
        name=Tasks.remove_orphaned_tasks, queue=settings.workflow_queue, max_retries=0
    )
    def remove_orphaned_tasks():
        lock = redlock(key=Tasks.remove_orphaned_tasks, auto_release_time=10.0)

        if not lock.acquire():
            return

        try:
            celery = get_celery_app()
            inspect = celery.control.inspect()
            tracker = OrphanedTasksTracker()
            task_ids = list(tracker.get_task_ids())

            if not task_ids:
                return

            active_workers = inspect.active() or {}
            workers = set(active_workers.keys())

            for task_id in task_ids:
                async_result = AsyncResult(id=task_id, app=celery)
                if async_result.state not in UNREADY_STATES:
                    continue
                if not async_result.info:
                    tracker.remove_task(task_id=task_id)
                    continue

                worker_name = async_result.info.get("hostname")
                if not worker_name:
                    tracker.remove_task(task_id)
                    continue

                if worker_name not in workers:
                    celery.backend.store_result(
                        task_id, Exception("Task worker is down"), state=FAILURE
                    )

        finally:
            if lock.locked():
                lock.release()
