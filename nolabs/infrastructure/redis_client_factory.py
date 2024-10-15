from contextlib import asynccontextmanager
from typing import Optional

import redis.asyncio as redis
from celery.bin.logtool import errors
from pottery import Redlock
from redis.client import Pipeline

from nolabs.infrastructure.log import logger
from nolabs.infrastructure.settings import settings

redis_client: Optional[redis.Redis] = None

class Redis:
    @classmethod
    @property
    def client(cls) -> redis.Redis:
        global redis_client

        if redis_client:
            return redis_client

        redis_client = redis.Redis.from_url(settings.celery_backend, decode_responses=True)
        return redis_client


def redlock(key: str, auto_release_time=100) -> Redlock:
    global redis_client

    if not redis_client:
        redis_client = Redis.client

    logger.info(f"Redlock captured", extra={"redlock_key": key, "redlock_auto_release_time": auto_release_time})
    return Redlock(key=key, masters={redis_client}, auto_release_time=auto_release_time)


@asynccontextmanager
async def use_redis_pipe():
    global redis_client

    if not redis_client:
        redis_client = Redis.client

    if isinstance(redis_client, Pipeline):
        yield
    else:
        original_client = redis_client
        pipe = redis_client.pipeline()
        redis_client = pipe
        try:
            yield
        finally:
            await pipe.execute()
            redis_client = original_client