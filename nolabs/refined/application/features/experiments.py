__all__ = [
    'GetExperimentsMetadataFeature',
    'CreateExperimentFeature'
]

import datetime
import uuid
from typing import List

from nolabs.refined.application.controllers.experiments.api_models import ExperimentMetadataResponse, \
    ChangeExperimentNameRequest
from nolabs.refined.domain.models.common import ExperimentId, ExperimentName
from nolabs.refined.domain.models.experiment import Experiment


def map_experiment_to_metadata(experiment: Experiment) -> ExperimentMetadataResponse:
    return ExperimentMetadataResponse(
        id=experiment.id.value,
        name=str(experiment.name),
        date=experiment.created_at
    )


class GetExperimentsMetadataFeature:
    def handle(self) -> List[ExperimentMetadataResponse]:
        experiments = Experiment.objects.all()

        result: List[ExperimentMetadataResponse] = []

        for experiment in experiments:
            result.append(map_experiment_to_metadata(experiment))

        return result


class CreateExperimentFeature:
    def handle(self) -> ExperimentMetadataResponse:
        experiment = Experiment(
            id=ExperimentId(uuid.uuid4()),
            name=ExperimentName('New experiment')
        )
        experiment.save()
        return map_experiment_to_metadata(experiment)


class DeleteExperimentFeature:
    def handle(self, experiment_id: uuid.UUID):
        assert experiment_id

        Experiment.objects.delete(id=ExperimentId(experiment_id))


class ChangeExperimentNameFeature:
    def handle(self, request: ChangeExperimentNameRequest):
        assert request

        experiment = Experiment.objects.get(id=ExperimentId(request.id))
        if experiment:
            experiment.set_name(request.name)
            experiment.save()
