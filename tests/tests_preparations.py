__all__ = [
    'mongo_disconnect',
    'mongo_connect'
]

import os

os.environ['NOLABS_ENVIRONMENT'] = 'unit'

import uuid
from datetime import datetime, timezone

from nolabs.domain.models.common import Experiment, ExperimentId, ExperimentName
from nolabs.infrastructure.settings import Settings

mongodb_connection = None


def mongo_connect():
    global mongodb_connection
    from nolabs.infrastructure.mongo_connector import mongo_connect
    settings = Settings.load()
    mongodb_connection = mongo_connect(settings.connection_string)



def mongo_disconnect():
    from nolabs.infrastructure.mongo_connector import mongo_disconnect  # type: ignore
    from mongoengine.connection import _get_db  # type: ignore

    settings = Settings.load()
    db_name = settings.connection_string.split('/')[-1]

    if mongodb_connection:
        mongodb_connection.drop_database(db_name)

    mongo_disconnect()


def seed_experiment() -> Experiment:
    experiment = Experiment(
        id=ExperimentId(uuid.uuid4()),
        name=ExperimentName('Test'),
        created_at=datetime.now(tz=timezone.utc)
    )
    experiment.save()
    return experiment
