from functools import lru_cache
from typing import Optional

import redis
from redis.client import Pipeline
from redis.lock import Lock

from nolabs.infrastructure.settings import settings


@lru_cache
def cached_client() -> redis.Redis:
    redis_client = redis.Redis.from_url(settings.celery_backend, decode_responses=True)
    return redis_client


class Redis:
    @classmethod
    @property
    def client(cls) -> redis.Redis:
        return cached_client()

    @classmethod
    async def disconnect(cls):
        client = cached_client()
        client.close()
        cached_client.cache_clear()


def redlock(key: str, blocking=False, auto_release_time=100) -> Optional[Lock]:
    return Lock(Redis.client, name=key, blocking=blocking, timeout=auto_release_time, blocking_timeout=auto_release_time)


def get_redis_pipe() -> Pipeline:
    return Redis.client.pipeline()
