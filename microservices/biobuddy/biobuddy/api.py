from fastapi import FastAPI, WebSocket, WebSocketDisconnect
from biobuddy.connection_manager import ConnectionManager
from biobuddy.services import send_message, send_message_async, invoke_action
from biobuddy.api_models import SendMessageToBioBuddyRequest, SendMessageToBioBuddyResponse, IsJobRunningResponse, \
    SendActionCallRequest
from typing import Dict
import json

app = FastAPI(
    title="Bio Buddy"
)

# Replace this with actual logic to manage stop tokens
stop_tokens: Dict[str, bool] = {}

connection_manager = ConnectionManager()

from biobuddy.job_state_manager import job_state_manager

@app.post("/send-message")
def predict(request: SendMessageToBioBuddyRequest) -> SendMessageToBioBuddyResponse:
    result = send_message(request)
    return result

@app.get("/job/{job_id}/is-running")
def is_job_running(job_id: str) -> IsJobRunningResponse:
    return IsJobRunningResponse(is_running=job_state_manager.is_job_running(job_id))

@app.get("/jobs/running")
def get_running_jobs():
    return {"running_jobs": job_state_manager.get_running_jobs()}


@app.post("/invoke-function")
async def invoke_function(
        request: SendActionCallRequest) -> SendMessageToBioBuddyResponse:
    action_text = request.action_text
    tools = request.tools
    return await invoke_action(action_text, tools)


@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    await connection_manager.connect(websocket)
    try:
        print("Original Websocket: ", websocket)
        while True:
            data = await websocket.receive_text()
            request_data = json.loads(data)

            if request_data.get('action') == 'stop':
                experiment_id = request_data.get('experiment_id')
                stop_tokens[experiment_id] = True
                await websocket.send_text(json.dumps({"reply_type": "final", "content": "<STOP>"}))
            else:
                request = SendMessageToBioBuddyRequest(**request_data)
                stop_tokens[request.experiment_id] = False
                async for message in send_message_async(request, stop_tokens, websocket):
                    if stop_tokens.get(request.experiment_id):
                        await websocket.send_text(json.dumps({"reply_type": "final", "content": "<STOP>"}))
                        break
                    await websocket.send_text(json.dumps(message.to_dict()))

    except WebSocketDisconnect:
        await connection_manager.disconnect(websocket)
