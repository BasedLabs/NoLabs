from celery import signals
from dotenv import load_dotenv

from nolabs.application import initialize
from nolabs.infrastructure.mongo_connector import mongo_connect
from nolabs.infrastructure.redis_client_factory import Redis

load_dotenv(".env")

from nolabs.infrastructure.log import initialize_logging, logger
from nolabs.infrastructure.settings import initialize_settings, settings
from nolabs.workflow.core.celery_tasks import register_workflow_celery_tasks
from nolabs.infrastructure.celery_app_factory import get_celery_app


@signals.worker_process_init.connect
def init_worker(**kwargs):
    mongo_connect()
    initialize()
    _ = Redis.client


def start():
    initialize_settings()
    initialize_logging()
    logger.info("Starting celery")
    app = get_celery_app()
    register_workflow_celery_tasks(app)
    app.worker_main(["worker", f"--concurrency={settings.celery_worker_concurrency}", "-P", settings.celery_worker_pool])
    return app

if __name__ == "__main__":
    start()