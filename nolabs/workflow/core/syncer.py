import asyncio
import time
import uuid

from nolabs.domain.exceptions import NoLabsException, ErrorCodes
from nolabs.infrastructure.celery_app_factory import get_celery_app
from workflow.core import Tasks


class Syncer:
    async def sync_graph(self, experiment_id: uuid.UUID, wait=False, timeout=604800):
        celery = get_celery_app()
        async_result = celery.send_task(name=Tasks.sync_graph_task, retry=False, kwargs={'experiment_id': experiment_id})
        if wait:
            start_time = time.time()
            while not async_result.ready():
                if time.time() - start_time > timeout:
                    raise NoLabsException(ErrorCodes.graph_scheduler_timeout)
                await asyncio.sleep(2.0)

