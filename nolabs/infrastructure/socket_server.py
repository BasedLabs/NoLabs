import asyncio
import queue
import threading
from typing import Any, Dict, Awaitable

import socketio

from nolabs.infrastructure.settings import settings

sio = socketio.AsyncServer(client_manager=socketio.AsyncRedisManager(settings.socketio_broker,
                                                           write_only=True))
sync_queue = queue.Queue()

def emit_event(name: str, room_id: str, data: Dict[str, Any]) -> Awaitable[None]:
    return sio.emit(event=name, room=room_id, data=data)

async def queue_worker():
    while True:
        if not sync_queue.empty():
            kwargs = sync_queue.get()
            await emit_event(**kwargs)
        else:
            await asyncio.sleep(0.1)

def start_async_worker():
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    loop.run_until_complete(queue_worker())

def enqueue_event(name: str, room_id: str, data: Dict[str, Any]):
    sync_queue.put({'name': name, 'room_id': room_id, 'data': data})

def start_queue_worker():
    threading.Thread(target=start_async_worker, daemon=True).start()