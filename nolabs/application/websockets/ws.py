__all__ = [
    'inject_websocket',
]

from typing import Dict, Any

from fastapi.websockets import WebSocket

websocket: WebSocket | None = None


def inject_websocket(ws: WebSocket):
    global websocket

    websocket = ws


async def send_websocket_message(event_type: str, message: Dict[str, Any]):
    global websocket

    if websocket:
        await websocket.send_json(data={
            'event_type': event_type,
            'message': message
        })
