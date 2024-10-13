from contextlib import asynccontextmanager
from typing import Optional

import redis.asyncio as redis
from pottery import Redlock

from nolabs.infrastructure.log import nolabs_logger as logger
from nolabs.infrastructure.settings import settings

redis_client: Optional[redis.Redis] = None

def get_redis_client() -> redis.Redis:
    global redis_client

    if redis_client:
        return redis_client

    redis_client = redis.Redis.from_url(settings().celery_backend, decode_responses=True)
    return redis_client


def redlock(key: str, auto_release_time=100) -> Redlock:
    global redis_client

    if not redis_client:
        redis_client = get_redis_client()

    logger.info(f"Redlock captured", extra={"redlock_key": key, "redlock_auto_release_time": auto_release_time})
    return Redlock(key=key, masters={redis_client}, auto_release_time=auto_release_time)

@asynccontextmanager
async def use_redis_pipe():
    global redis_client

    if not redis_client:
        redis_client = get_redis_client()

    original_client = redis_client
    pipe = redis_client.pipeline()
    redis_client = pipe
    try:
        yield
    finally:
        await pipe.execute()
        redis_client = original_client