import redis
from pottery import Redlock

from nolabs.infrastructure.settings import settings

host = settings.celery_backend.split('/')[-2].split(':')[0]
port = int(settings.celery_backend.split('/')[-2].split(':')[1])

redis_client = redis.StrictRedis(host=host, port=port, db=1)


def acquire_redlock(key: str, timeout=5):
    return Redlock(key=key, masters={redis_client},
                   raise_on_redis_errors=True,
                   context_manager_blocking=True,
                   context_manager_timeout=timeout)
