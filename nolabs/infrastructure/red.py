from contextlib import asynccontextmanager

import redis.asyncio as redis
from pottery import Redlock

from nolabs.infrastructure.log import logger
from nolabs.infrastructure.settings import settings

redis_client = redis.Redis.from_url(settings.celery_backend)


def redlock(key: str, auto_release_time=100) -> Redlock:
    logger.info(f"Redlock captured", extra={"redlock_key": key, "redlock_auto_release_time": auto_release_time})
    return Redlock(key=key, masters={redis_client}, auto_release_time=auto_release_time)

@asynccontextmanager
async def use_redis_pipe():
    global redis_client
    original_client = redis_client
    pipe = redis_client.pipeline()
    redis_client = pipe
    try:
        yield
    finally:
        await pipe.execute()
        redis_client = original_client