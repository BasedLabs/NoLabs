from nolabs.domain.experiment import ExperimentId
from nolabs.api_models.gene_ontology import GetExperimentResponse, ExperimentMetadataResponse
from nolabs.features.gene_ontology.services.file_management import FileManagement


class GetExperimentFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    def handle(self, id: str) -> GetExperimentResponse:
        assert id

        experiment_id = ExperimentId(id)
        metadata = self._file_management.get_experiment_metadata()
        data = self._file_management.get_experiment_data(experiment_id)
        return GetExperimentResponse(
            metaData=ExperimentMetadataResponse(
                id=metadata.id.value,
                name=metadata.name.value,
                date=metadata.date
            ),
            data=data
        )

