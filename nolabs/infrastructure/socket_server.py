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

    async def queue_worker(self):
        while True:
            if not self.sync_queue.empty():
                kwargs = self.sync_queue.get()
                await self.emit_event(**kwargs)
            else:
                await asyncio.sleep(0.1)

    def start_async_worker(self):
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        loop.run_until_complete(self.queue_worker())

    def enqueue_event(self, name: str, room_id: str, data: Dict[str, Any]):
        self.sync_queue.put({"name": name, "room_id": room_id, "data": data})

    def start_queue_worker(self):
        asyncio.create_task(self.start_async_worker())

    def disconnect(self):
        self.client.redis.close()
        get_socket_server.cache_clear()


@lru_cache
def get_socket_server():
    return SocketServer()
