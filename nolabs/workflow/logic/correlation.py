import uuid
from typing import List, Dict, Any

from nolabs.infrastructure.red import redis_client


def _make_component_cid(id: str | uuid.UUID) -> str:
    return f"correlation_id:component:{str(id)}"


def _make_job_cid(id: uuid.UUID) -> str:
    return f"correlation_id:job:{str(id)}"


def unpack_cid(cid: str) -> uuid.UUID:
    return uuid.UUID(cid.split(':')[-1])


async def _assign_correlation_id(correlation_id: str, step: str) -> str:
    """
    Call this RIGHT BEFORE next task in chain, return immediately after
    :return: celery_task_id
    """
    await redis_client.delete(correlation_id)
    celery_task_id = str(uuid.uuid4())
    await redis_client.rpush(correlation_id, celery_task_id)
    return celery_task_id


async def _assign_correlation_id(correlation_id: str, step: str) -> str:
    """
    Call this RIGHT BEFORE next task in chain, return immediately after
    :return: celery_task_id
    """
    await redis_client.delete(correlation_id)
    celery_task_id = str(uuid.uuid4())
    await redis_client.hset(correlation_id, step, celery_task_id)
    return celery_task_id


async def _get_correlation_data(correlation_id: str) -> Dict[str, Any]:
    d = await redis_client.hgetall(name=correlation_id)
    return {k.decode('utf-8'): v.decode('utf-8') if isinstance(v, bytes) else v for k, v in d.items()}


async def _clear_correlation(cid: str):
    await redis_client.delete(cid)


async def _all_components_cids() -> List[str]:
    _, keys = await redis_client.scan(match="correlation_id:component:*")
    return keys


async def _all_jobs_cids() -> List[str]:
    _, keys = await redis_client.scan(match="correlation_id:job:*")
    return keys


async def _all_celery_tasks(cid: str) -> List[str]:
    return await redis_client.lrange(cid, 0, -1)
