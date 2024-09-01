import logging
from settings import settings

mode = settings.mode

logger = logging.getLogger(__name__)


if mode == 'celery':
    from worker import app
    app.worker_main([
        'worker',
        f'--loglevel={settings.application_loglevel}'
    ])

if mode == 'fastapi':
    import uvicorn

    from fastapi_api import app

    uvicorn.run(app, host=settings.fastapi_host, port=settings.fastapi_port)
