__all__ = [
    'GetExperimentsMetadataFeature',
    'CreateExperimentFeature'
]

import uuid
from typing import List

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.application.use_cases.experiments.api_models import ExperimentMetadataResponse, \
    UpdateExperimentRequest
from nolabs.domain.models.common import ExperimentId, ExperimentName, Experiment


def map_experiment_to_metadata(experiment: Experiment) -> ExperimentMetadataResponse:
    return ExperimentMetadataResponse(
        id=experiment.id,
        name=str(experiment.name),
        date=experiment.created_at
    )


class GetExperimentsMetadataFeature:
    def handle(self) -> List[ExperimentMetadataResponse]:
        experiments: List[Experiment] = Experiment.objects.all()

        result: List[ExperimentMetadataResponse] = []

        for experiment in experiments:
            result.append(map_experiment_to_metadata(experiment))

        return result


class CreateExperimentFeature:
    def handle(self) -> ExperimentMetadataResponse:
        experiment = Experiment.create(
            id=ExperimentId(uuid.uuid4()),
            name=ExperimentName('New experiment')
        )
        experiment.save()
        return map_experiment_to_metadata(experiment)


class DeleteExperimentFeature:
    async def handle(self, experiment_id: uuid.UUID):
        assert experiment_id

        experiment: Experiment = Experiment.objects.with_id(experiment_id)
        if experiment:
            await experiment.delete()


class UpdateExperimentFeature:
    def handle(self, request: UpdateExperimentRequest):
        assert request

        experiment = Experiment.objects.with_id(request.id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        if request.name:
            experiment.set_name(ExperimentName(request.name))

        experiment.save()
