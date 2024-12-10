from typing import Any, List

from dotenv import load_dotenv
from starlette.responses import HTMLResponse, FileResponse, StreamingResponse

load_dotenv(".env")

from settings import settings
import application_paper_qa as application
from api_models import ChatRequest, ChatResponse, SetOpenAiApiKey
from fastapi import FastAPI
import uvicorn

app = FastAPI(title="Arxiv abstracts AI")

async def chat_generator(request: ChatRequest):
    response: ChatResponse
    async for state in application.inference(request=request):
        yield state.model_dump_json().encode("utf-8") + b"\n"

@app.post("/api/chat")
async def inference(request: ChatRequest) -> StreamingResponse:
    return StreamingResponse(chat_generator(request=request), media_type="application/json")

@app.get('/api/chat/history')
async def history() -> List[ChatResponse]:
    return await application.history()

@app.delete('/api/chat/history')
async def clear_history():
    application.clear_history()

@app.post('/api/openai-api-key')
async def set_api_key(request: SetOpenAiApiKey):
    application.set_openai_api_key(api_key=request.api_key)

@app.get("/chat", response_class=HTMLResponse)
async def serve_chat_page() -> Any:
    return FileResponse("simple-ui.html")

if __name__ == "__main__":
    uvicorn.run(app, host=settings.fastapi_host, port=settings.fastapi_port)
