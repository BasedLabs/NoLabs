import datetime

from celery import Celery
from celery.result import AsyncResult
from celery.states import FAILURE, READY_STATES

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
            tasks = list(OrphanedTasksTracker.get_task_data())

            if not tasks:
                return

            active_workers = inspect.active() or {}
            active_queues = inspect.active_queues()
            workers = set(active_workers.keys())
            queues = set()
            for l in active_queues.values():
                for queue in l:
                    queues.add(queue['name'])

            for task in tasks:
                task_id = task['task_id']
                task_queue = task['queue']
                timestamp = task['timestamp']
                async_result = AsyncResult(id=task_id, app=celery)

                if async_result.state in READY_STATES:
                    OrphanedTasksTracker.remove_task(task_id)
                    continue

                if datetime.datetime.now(datetime.UTC).timestamp() - timestamp < 60.0:
                    continue

                if async_result.info:
                    worker_name = async_result.info.get('hostname')
                    if worker_name and worker_name not in workers:
                        celery.backend.store_result(
                            task_id, Exception("Task worker is down (failed in progress)"), state=FAILURE
                        )
                        OrphanedTasksTracker.remove_task(task_id)
                        continue

                if task_queue not in queues:
                    celery.backend.store_result(
                        task_id, Exception("Task worker is down (not started)"), state=FAILURE
                    )
                    OrphanedTasksTracker.remove_task(task_id)
                    continue
        finally:
            if lock.owned():
                lock.release()
