from dotenv import load_dotenv
from log import logger
from settings import settings

mode = settings.mode

load_dotenv(".env")

if mode == "celery":
    from worker import app

    app.worker_main(["worker", f"--concurrency={settings.celery_worker_concurrency}"])

if mode == "fastapi":
    import uvicorn
    from fastapi_api import app

    uvicorn.run(app, host=settings.fastapi_host, port=settings.fastapi_port)
