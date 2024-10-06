import uuid
from typing import List

from nolabs.infrastructure.red import redis_client


def _make_correlation_id(id: str | uuid.UUID) -> str:
    return f"correlation_id:{str(id)}"

def _unpack_correlation_id(correlation_id: str) -> uuid.UUID:
    return uuid.UUID(correlation_id.split(':')[-1])

async def _assign_correlation_id(correlation_id: str) -> str:
    """
    :return: celery_task_id
    """
    await redis_client.delete(correlation_id)
    celery_task_id = str(uuid.uuid4())
    await redis_client.rpush(correlation_id, celery_task_id)
    return celery_task_id

async def _all_correlation_ids() -> List[str]:
    _, keys = await redis_client.scan(match="correlation_id:*")
    return keys

async def _all_celery_tasks(cid: str) -> List[str]:
    return await redis_client.lrange(cid, 0, -1)