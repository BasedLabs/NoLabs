from dotenv import load_dotenv
from log import logger
from settings import settings

load_dotenv(".env")

import uvicorn
from fastapi_api import app

logger.info("Starting fastapi")
uvicorn.run(app, host=settings.fastapi_host, port=settings.fastapi_port)
