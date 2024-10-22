from typing import Optional

import redis.asyncio as redis
from pottery import Redlock

from nolabs.infrastructure.log import logger
from nolabs.infrastructure.settings import settings

class Redis:
    @classmethod
    @property
    def client(cls) -> redis.Redis:
        redis_client = redis.Redis.from_url(settings.celery_backend, decode_responses=True)
        return redis_client


def redlock(key: str, auto_release_time=100) -> Redlock:
    redis_client = Redis.client
    logger.info(f"Redlock captured", extra={"redlock_key": key, "redlock_auto_release_time": auto_release_time})
    return Redlock(key=key, masters={redis_client}, auto_release_time=auto_release_time)


def get_redis_pipe():
    return Redis.client.pipeline()
