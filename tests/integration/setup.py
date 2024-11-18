import multiprocessing
import os
import time
import tracemalloc
import unittest
import uuid
from pathlib import Path
from typing import Optional

from testcontainers.mongodb import MongoDbContainer
from testcontainers.redis import RedisContainer

from nolabs.infrastructure.celery_app_factory import get_celery_app
from nolabs.infrastructure.mongo_connector import mongo_connect, mongo_disconnect
from nolabs.infrastructure.redis_client_factory import cached_client
from nolabs.infrastructure.settings import initialize_settings
from nolabs.workflow.core.celery_tasks import register_workflow_celery_tasks
from nolabs.workflow.core.component import ComponentTypeFactory
from nolabs.workflow.worker import start

redis_container: Optional[RedisContainer] = None
mongo_container: Optional[MongoDbContainer] = None

tracemalloc.start()


class GlobalSetup(unittest.IsolatedAsyncioTestCase):
    worker_thread = None
    setup_called = False

    @classmethod
    def setUpClass(cls):
        """Start Redis and MongoDB containers before any test class runs."""
        global redis_container, mongo_container

        # Start Redis container
        redis_container = RedisContainer()
        redis_container.with_exposed_ports(6379)
        redis_container.with_bind_ports(6379, 22555)
        redis_container.start()

        # Start MongoDB container
        mongo_container = MongoDbContainer(
            username="admin", password="admin", dbname="nolabs"
        )
        mongo_container.with_exposed_ports(27017)
        mongo_container.with_bind_ports(27017, 22556)
        mongo_container.start()

        # Set environment variables
        redis_host = redis_container.get_container_host_ip()
        redis_port = redis_container.get_exposed_port(6379)
        mongo_host = mongo_container.get_container_host_ip()
        mongo_port = mongo_container.get_exposed_port(27017)

        os.environ["NOLABS_BIOBUDDY_HOST"] = "http://mock-biobuddy-host"
        os.environ["NOLABS_EXTERNAL_QUERY_HOST"] = "http://mock-query-host"
        os.environ["NOLABS_ESMFOLD_LIGHT_HOST"] = "http://mock-esmfold-host"
        os.environ["NOLABS_REINVENT_HOST"] = "http://mock-reinvent-host"
        os.environ["NOLABS_DIFFDOCK_HOST"] = "http://mock-diffdock-host"
        os.environ["NOLABS_CELERY_BROKER"] = f"redis://{redis_host}:{redis_port}/0"
        os.environ["NOLABS_CELERY_BACKEND"] = f"redis://{redis_host}:{redis_port}/0"
        os.environ["NOLABS_CELERY_WORKER_STATE_DB"] = f"/tmp/{str(uuid.uuid4())}.db"
        os.environ["NOLABS_CELERY_EAGER"] = "False"
        os.environ["NOLABS_CONNECTION_STRING"] = (
            f"mongodb://admin:admin@{mongo_host}:{mongo_port}/nolabs?authSource=admin"
        )
        os.environ["NOLABS_SOCKETIO_BROKER"] = f"redis://{redis_host}:{redis_port}/0"
        os.environ["NOLABS_CELERY_WORKER_POOL"] = "prefork"
        os.environ["NOLABS_REINVENT_DIRECTORY"] = str(Path("/tmp/reinvent/dir"))
        os.environ["NOLABS_BLAST_EMAIL"] = "mock-email@example.com"
        os.environ["NOLABS_WORKFLOW_VERSION"] = "2"
        os.environ["NOLABS_MODE"] = "united"
        os.environ["NOLABS_ENVIRONMENT"] = "test"
        os.environ["NOLABS_LOGGING_LEVEL"] = "ERROR"

        initialize_settings()
        # initialize_logging()

        register_workflow_celery_tasks(get_celery_app())

        time.sleep(5.0)

    @classmethod
    def tearDownClass(cls):
        global redis_container, mongo_container

        if redis_container:
            redis_container.stop()
        if mongo_container:
            mongo_container.stop()

    def setUp(self):
        mongo_connect()

    def tearDown(self):
        ComponentTypeFactory.clear()

        app = get_celery_app()
        app.control.shutdown()
        get_celery_app.cache_clear()
        cached_client.cache_clear()

        mongo_disconnect()

        if self.worker_process:
            self.worker_process.terminate()
            self.worker_process = None

    def spin_up_celery(self):
        # self.worker_thread = threading.Thread(target=start)
        # self.worker_thread.daemon = True
        # self.worker_thread.start()

        self.worker_process = multiprocessing.Process(target=start)
        self.worker_process.daemon = (
            True  # If needed, the process will exit when the main program does
        )
        self.worker_process.start()
