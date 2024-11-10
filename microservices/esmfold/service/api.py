from dotenv import load_dotenv

load_dotenv(".env")

from log import logger
from settings import settings




import uvicorn
from fastapi_api import app

logger.info("Starting fastapi")
uvicorn.run(app, host=settings.fastapi_host, port=settings.fastapi_port)
