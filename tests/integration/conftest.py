import asyncio
import os
import threading
from pathlib import Path
from typing import Optional

import pytest
import pytest_asyncio
from testcontainers.mongodb import MongoDbContainer
from testcontainers.redis import RedisContainer

# Global variables to hold container references
redis_container: Optional[RedisContainer] = None
mongo_container: Optional[MongoDbContainer] = None


#@pytest_asyncio.fixture(scope='function')
#def event_loop(request):
#    """Create an instance of the default event loop for each test case."""
#    loop = asyncio.get_event_loop_policy().new_event_loop()
#    yield loop
#    loop.close()


def start_containers():
    """Start Redis and MongoDB containers globally."""
    global redis_container, mongo_container

    # Start Redis container
    redis_container = RedisContainer()
    redis_container.start()

    # Start MongoDB container
    mongo_container = MongoDbContainer(username="admin", password="admin", dbname="nolabs")
    mongo_container.start()


def stop_containers():
    """Stop Redis and MongoDB containers."""
    global redis_container, mongo_container
    if redis_container:
        redis_container.stop()
    if mongo_container:
        mongo_container.stop()


@pytest.hookimpl(tryfirst=True)
def pytest_configure(config):
    """
    This hook is run before any test modules are imported.
    It's perfect for setting up global environment variables and containers.
    """
    start_containers()  # Start containers before tests are imported

    redis_host = redis_container.get_container_host_ip()
    redis_port = redis_container.get_exposed_port(6379)
    mongo_host = mongo_container.get_container_host_ip()
    mongo_port = mongo_container.get_exposed_port(27017)

    # Set environment variables before any test modules are imported
    os.environ["NOLABS_BIOBUDDY_HOST"] = "http://mock-biobuddy-host"
    os.environ["NOLABS_EXTERNAL_QUERY_HOST"] = "http://mock-query-host"
    os.environ["NOLABS_ESMFOLD_LIGHT_HOST"] = "http://mock-esmfold-host"
    os.environ["NOLABS_REINVENT_HOST"] = "http://mock-reinvent-host"
    os.environ["NOLABS_DIFFDOCK_HOST"] = "http://mock-diffdock-host"
    os.environ["NOLABS_CELERY_BROKER"] = f"redis://{redis_host}:{redis_port}/0"
    os.environ["NOLABS_CELERY_BACKEND"] = f"redis://{redis_host}:{redis_port}/0"
    os.environ["NOLABS_CELERY_EAGER"] = "False"
    os.environ["NOLABS_CONNECTION_STRING"] = f"mongodb://admin:admin@{mongo_host}:{mongo_port}/nolabs?authSource=admin"
    os.environ["NOLABS_SOCKETIO_BROKER"] = f"redis://{redis_host}:{redis_port}/0"
    os.environ["NOLABS_CELERY_WORKER_POOL"] = "threads"
    os.environ["NOLABS_REINVENT_DIRECTORY"] = str(Path("/tmp/reinvent/dir"))
    os.environ["NOLABS_BLAST_EMAIL"] = "mock-email@example.com"
    os.environ["NOLABS_WORKFLOW_VERSION"] = "2"
    os.environ["NOLABS_MODE"] = "united"
    os.environ["NOLABS_ENVIRONMENT"] = "test"
    os.environ["NOLABS_ENABLE_STRUCTURED_LOGGING"] = "False"
    os.environ["NOLABS_LOGGING_LEVEL"] = "ERROR"


def pytest_unconfigure(config):
    """
    This hook is called at the end of the test session.
    It's used to clean up any global resources, such as stopping containers.
    """
    stop_containers()


@pytest_asyncio.fixture
def prefork_celery_worker():
    """Fixture to run Celery worker in a separate process for testing."""
    from workflow.worker import start
#
    #worker_process = multiprocessing.Process(target=start)
    #worker_process.start()
    ## time.sleep(5)  # Allow some time for the worker to start
#
    #yield  # Tests will run while the worker is active
#
    ## Clean up the worker after tests
    #worker_process.terminate()
    #worker_process.join()
    worker_thread = threading.Thread(target=start)
    worker_thread.daemon = True

    # Start the worker thread
    worker_thread.start()

    # Optionally, wait a few seconds to ensure the worker has started

    yield  # This is where your tests will run while the worker is active

    # Clean up the worker thread after tests
    # You should ensure that your worker has a way to stop gracefully
    # as threads do not have `terminate()` like processes
    # This assumes that `start()` handles stopping based on some condition
    # (e.g., using a signal or flag for graceful shutdown).

    # Ensure the worker thread finishes
    #worker_thread.join()

@pytest_asyncio.fixture(autouse=True)
def clear_type_factory():
    from nolabs.workflow.core.component import ComponentTypeFactory

    """Fixture to start asserts before and after a test is run"""
    yield

    ComponentTypeFactory.clear()