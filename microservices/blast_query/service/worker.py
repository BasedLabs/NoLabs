from typing import Dict, Any

from dotenv import load_dotenv

load_dotenv(".env")

import uuid
from celery import Celery
from settings import settings
from api_models import BlastType
import redis

app = Celery(
    __name__, backend=settings.redis_url, broker=settings.redis_url
)
app.conf.update(
    task_default_queue=settings.celery_worker_queue,
    task_track_started=True,
    task_acks_late=True,
    task_reject_on_worker_lost=True,
    worker_state_db="/opt/celery-state.db",
    accept_content=["application/json", "application/x-python-serialize"],
    broker_transport_options={"heartbeat": 10},
    broker_connection_retry_on_startup=True,
    worker_send_task_events=True
)

def ensure_redis_available():
    try:
        client = redis.from_url(settings.redis_url, socket_connect_timeout=5)
        if client.ping():
            return
    except redis.ConnectionError:
        raise Exception('Redis is not connected. Fix REDIS_URL in corresponding .env file.')

@app.task(time_limit=720, name="run_blast", rate_limit="5/m")
def run_blast(param: Dict[str, Any]) -> Dict[str, Any]:
    import application

    result = application.run_blast(
        program=BlastType.blastp,
        query=param['query'],
        descriptions=param['descriptions'],
        alignments=param['alignments'],
        hitlist_size=param['hitlist_size'],
        expect=param['expect']
    )

    return result

if __name__ == "__main__":
    ensure_redis_available()
    app.worker_main(
        [
            "worker",
            f"--concurrency={settings.celery_worker_concurrency}",
            "-E",
            "-n",
            f"{settings.celery_worker_queue}-{str(uuid.uuid4())}",
            "--loglevel", settings.logging_level
        ]
    )
