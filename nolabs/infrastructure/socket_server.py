import asyncio
import queue
from functools import lru_cache
from typing import Any, Awaitable, Dict

import socketio

from nolabs.infrastructure.settings import settings


class SocketServer:
    def __init__(self):
        self.client = socketio.RedisManager(settings.socketio_broker, write_only=True)
        self.sio = socketio.Server(client_manager=self.client)
        self.sync_queue = queue.Queue()

    def emit_event(
        self, name: str, room_id: str, data: Dict[str, Any]
    ) -> Awaitable[None]:
        return self.sio.emit(event=name, room=room_id, data=data)

    def disconnect(self):
        self.client.redis.close()
        get_socket_server.cache_clear()


@lru_cache
def get_socket_server():
    return SocketServer()
