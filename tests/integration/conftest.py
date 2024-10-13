import multiprocessing
import time
from pathlib import Path

import pytest
from testcontainers.mongodb import MongoDbContainer
from testcontainers.redis import RedisContainer

from infrastructure.log import init_logging
from infrastructure.settings import init_settings


@pytest.fixture
def redis_container():
    """Set up a Redis container."""
    with RedisContainer() as redis:
        yield redis


@pytest.fixture
def mongo_container():
    """Set up a MongoDB container."""
    with MongoDbContainer() as mongo:
        yield mongo


@pytest.fixture
def setup(monkeypatch, redis_container, mongo_container):
    redis_host = redis_container.get_container_host_ip()
    redis_port = redis_container.get_exposed_port(6379)

    mongo_host = mongo_container.get_container_host_ip()
    mongo_port = mongo_container.get_exposed_port(27017)

    monkeypatch.setenv("NOLABS_BIOBUDDY_HOST", "http://mock-biobuddy-host")
    monkeypatch.setenv("NOLABS_EXTERNAL_QUERY_HOST", "http://mock-query-host")
    monkeypatch.setenv("NOLABS_ESMFOLD_LIGHT_HOST", "http://mock-esmfold-host")
    monkeypatch.setenv("NOLABS_REINVENT_HOST", "http://mock-reinvent-host")
    monkeypatch.setenv("NOLABS_DIFFDOCK_HOST", "http://mock-diffdock-host")
    monkeypatch.setenv("NOLABS_CELERY_BROKER", f"redis://{redis_host}:{redis_port}/0")
    monkeypatch.setenv("NOLABS_CELERY_BACKEND", f"redis://{redis_host}:{redis_port}/0")
    monkeypatch.setenv("NOLABS_CELERY_EAGER", "False")
    monkeypatch.setenv("NOLABS_CONNECTION_STRING", f"mongodb://{mongo_host}:{mongo_port}/nolabs")
    monkeypatch.setenv("NOLABS_SOCKETIO_BROKER", f"redis://{redis_host}:{redis_port}/0")
    monkeypatch.setenv("NOLABS_REINVENT_DIRECTORY", str(Path("/tmp/reinvent/dir")))
    monkeypatch.setenv("NOLABS_BLAST_EMAIL", "mock-email@example.com")
    monkeypatch.setenv("NOLABS_WORKFLOW_VERSION", "2")
    monkeypatch.setenv("NOLABS_MODE", "united")
    monkeypatch.setenv("NOLABS_ENVIRONMENT", "test")
    monkeypatch.setenv("NOLABS_ENABLE_STRUCTURED_LOGGING", "False")
    monkeypatch.setenv("NOLABS_LOGGING_LEVEL", "ERROR")

    init_settings()
    init_logging()

@pytest.fixture
def prefork_celery_worker():
    """Fixture to run Celery worker in a separate process for testing."""
    from workflow.worker import start

    worker_process = multiprocessing.Process(target=start)
    worker_process.start()
    time.sleep(5)  # Allow some time for the worker to start

    yield  # Tests will run while the worker is active

    # Clean up the worker after tests
    worker_process.terminate()
    worker_process.join()