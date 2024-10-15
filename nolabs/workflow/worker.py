from dotenv import load_dotenv


load_dotenv(".env")

from infrastructure.log import worker_logger as logger
from infrastructure.celery_app_factory import get_celery_app
from infrastructure.settings import settings
import workflow.core.celery_tasks # noqa: F401 Used for initializing celery tasks


def start():
    logger.info("Starting celery")
    app = get_celery_app()
    app.worker_main(["worker", f"--concurrency={settings.celery_worker_concurrency}"])

if __name__ == "__main__":
    start()