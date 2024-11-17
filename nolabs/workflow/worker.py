from dotenv import load_dotenv

load_dotenv(".env")

import uuid

from asgiref.sync import async_to_sync
from celery import signals

from nolabs.application import initialize
from nolabs.infrastructure.celery_app_factory import get_celery_app
from nolabs.infrastructure.mongo_connector import mongo_connect, mongo_disconnect
from nolabs.infrastructure.redis_client_factory import rd
from nolabs.infrastructure.settings import initialize_settings, settings
from nolabs.infrastructure.socket_server import get_socket_server
from nolabs.workflow.core import Tasks
from nolabs.workflow.core.celery_tasks import register_workflow_celery_tasks
from nolabs.workflow.monitoring.celery_tasks import register_monitoring_celery_tasks


@signals.worker_ready.connect
def task_prerun(**kwargs):
    initialize()
    mongo_connect()


@signals.worker_shutting_down.connect
def task_postrun(**kwargs):
    async def _():
        rd.close()
        mongo_disconnect()
        get_socket_server().disconnect()

    async_to_sync(_)()


def start(beat=False):
    initialize_settings()
    app = get_celery_app()
    app.conf.update(
        beat_schedule={
            Tasks.sync_graphs_task: {
                "task": Tasks.sync_graphs_task,
                "schedule": 1.0,
                "args": (),
            },
            Tasks.remove_orphaned_tasks: {
                "task": Tasks.remove_orphaned_tasks,
                "schedule": 30.0,
                "args": (),
            },
        },
        task_acks_late=True,
        broker_transport_options={"heartbeat": 10},
        task_soft_time_limit=10,
        task_time_limit=15,
        redbeat_lock_timeout=20,
        broker_connection_retry_on_startup=True,
        task_track_started=True,
        task_reject_on_worker_lost=True,
        accept_content=["application/json", "application/x-python-serialize"],
        worker_send_task_events=True
    )
    app.autodiscover_tasks(force=True)
    register_workflow_celery_tasks(app)
    register_monitoring_celery_tasks(app)
    args = [
        "worker",
        f"--concurrency={settings.celery_worker_concurrency}",
        "-P",
        settings.celery_worker_pool,
        "-n",
        f"workflow-{str(uuid.uuid4())}",
        "--scheduler",
        "redbeat.RedBeatScheduler"
    ]
    if beat:
        args.append("-B")
    app.worker_main(args)
    return app


if __name__ == "__main__":
    start(beat=True)
