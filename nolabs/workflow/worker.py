from dotenv import load_dotenv


load_dotenv(".env")

from nolabs.infrastructure.log import initialize_logging, logger
from nolabs.infrastructure.settings import initialize_settings, settings
from nolabs.workflow.core.celery_tasks import register_workflow_celery_tasks
from nolabs.infrastructure.celery_app_factory import get_celery_app


def start():
    initialize_settings()
    initialize_logging()
    logger.info("Starting celery")
    app = get_celery_app()
    register_workflow_celery_tasks(app)
    app.worker_main(["worker", f"--concurrency={settings.celery_worker_concurrency}", "-P", "threads"])
    return app

if __name__ == "__main__":
    start()