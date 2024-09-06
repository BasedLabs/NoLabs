import json
from typing import Optional, Any, Dict

import redis

from nolabs.infrastructure.settings import settings


class WebsocketsQueue:
    def __init__(self, redis_host: str, redis_port: int):
        self.redis_client = redis.Redis(host=redis_host, port=redis_port, db=0)

    def read_last(self) -> Optional[Dict[str, Any]]:
        item = self.redis_client.rpop('websockets')
        if item:
            return json.loads(item)
        return None

    def write(self, data: Dict[str, Any]):
        self.redis_client.lpush('websockets', json.dumps(data))

    def clear_db(self):
        self.redis_client.flushdb()


websockets_queue = WebsocketsQueue(settings.redis_host, settings.redis_port)
