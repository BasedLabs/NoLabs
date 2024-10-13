from typing import Optional

from celery import Celery

from nolabs.infrastructure.settings import settings

app: Optional[Celery] = None


def get_celery_app() -> Celery:
    global app

    if app:
        return app

    app = Celery(
        __name__, broker=settings.celery_broker, backend=settings.celery_backend
    )
    app.conf.update(
        enable_utc=True,
        task_default_queue='workflow',
        task_track_started=True,
        worker_send_task_events=True,
        task_acks_late=True,
        task_reject_on_worker_lost=True,
        worker_state_db="/tmp/celery-state.db",
        accept_content=["application/json", "application/x-python-serialize"],
        task_always_eager=False,
        task_eager_propagates=False
    )
    app.autodiscover_tasks(force=True)

    return app