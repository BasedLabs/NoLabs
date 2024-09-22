from typing import Any, Dict

import socketio

from nolabs.infrastructure.settings import settings

sio = socketio.Server(client_manager=socketio.RedisManager(settings.socketio_broker))


def emit_event(name: str, id: str, data: Dict[str, Any]):
    sio.emit(event=name, room=id, data=data)
