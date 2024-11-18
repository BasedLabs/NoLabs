from functools import lru_cache
from typing import Optional

import redis
from redis.client import Pipeline
from redis.lock import Lock

from nolabs.infrastructure.settings import settings


@lru_cache
def cached_client() -> redis.Redis:
    redis_client = redis.Redis.from_url(settings.redis_url, decode_responses=True)
    return redis_client


class RedisProxy:
    def __getattr__(self, name):
        return getattr(cached_client(), name)

    def __call__(self):
        return cached_client()


rd: redis.Redis = RedisProxy()  # type: ignore


def redlock(key: str, blocking=False, auto_release_time=100) -> Optional[Lock]:
    return Lock(
        rd,
        name=key,
        blocking=blocking,
        timeout=auto_release_time,
        blocking_timeout=auto_release_time,
    )


def get_redis_pipe() -> Pipeline:
    return rd.pipeline()

def ensure_redis_available():
    try:
        client = redis.from_url(settings.redis_url, socket_connect_timeout=5)
        if client.ping():
            return
    except redis.ConnectionError as e:
        raise Exception('Redis is not connected. Fix REDIS_URL in corresponding .env file.')

