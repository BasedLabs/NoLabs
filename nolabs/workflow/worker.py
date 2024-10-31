from asgiref.sync import async_to_sync
from celery import signals, Celery
from celery.schedules import crontab
from dotenv import load_dotenv

from nolabs.application import initialize
from nolabs.infrastructure.mongo_connector import mongo_connect, mongo_disconnect
from nolabs.infrastructure.redis_client_factory import Redis
from nolabs.infrastructure.socket_server import get_socket_server
from nolabs.workflow.core import Tasks

load_dotenv(".env")

from nolabs.infrastructure.log import initialize_logging, logger
from nolabs.infrastructure.settings import initialize_settings, settings
from nolabs.workflow.core.celery_tasks import register_workflow_celery_tasks
from nolabs.infrastructure.celery_app_factory import get_celery_app


@signals.task_prerun.connect
def task_prerun(**kwargs):
    initialize()
    mongo_connect()


@signals.task_postrun.connect
def task_postrun(**kwargs):
    async def _():
        await Redis.disconnect()
        Redis.clear_cache()
        mongo_disconnect()
        get_socket_server().disconnect()

    async_to_sync(_)()


def start(beat=False):
    initialize_settings()
    initialize_logging()
    logger.info("Starting celery")
    app = get_celery_app()
    app.conf.update(
        beat_schedule={
            Tasks.sync_graphs_task: {
                'task': Tasks.sync_graphs_task,
                'schedule': 1.0,
                'args': ()
            }
        }
    )
    app.autodiscover_tasks(force=True)
    register_workflow_celery_tasks(app)
    args = ["worker", f"--concurrency={settings.celery_worker_concurrency}", "-P", settings.celery_worker_pool,
         f"--loglevel={settings.logging_level}"]
    if beat:
        args.append('-B')
    app.worker_main(args)
    return app


if __name__ == "__main__":
    start(beat=True)
