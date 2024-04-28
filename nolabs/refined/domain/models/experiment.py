import datetime

from mongoengine import Document, DateTimeField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import ExperimentId, ExperimentName
from nolabs.refined.infrastructure.mongo_fields import ValueObjectUUIDField, ValueObjectStringField
from nolabs.seedwork.domain.entities import Entity


class Experiment(Document, Entity):
    id: ExperimentId = ValueObjectUUIDField(db_field='_id', primary_key=True, required=True)
    name: ExperimentName = ValueObjectStringField(required=True)
    created_at: datetime.datetime = DateTimeField()

    def __init__(self, id: ExperimentId, name: ExperimentName, *args, **values):
        if not id:
            raise NoLabsException('Experiment id is empty', ErrorCodes.invalid_experiment_id)

        if not name:
            raise NoLabsException('Experiment name is empty', ErrorCodes.invalid_experiment_name)

        created_at = datetime.datetime.now(tz=datetime.timezone.utc)

        super().__init__(id=id, name=name, created_at=created_at, *args, **values)

    def set_name(self, name: ExperimentName):
        if not name:
            raise NoLabsException('Name cannot be empty', ErrorCodes.invalid_experiment_name)

        self.name = name

    def save(self, *args, **kwargs):
        super().save(*args, **kwargs)
