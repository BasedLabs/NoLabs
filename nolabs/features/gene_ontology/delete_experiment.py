from nolabs.domain.experiment import ExperimentId
from nolabs.features.gene_ontology.services.file_management import FileManagement


class DeleteExperimentFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    def handle(self, id: str):
        assert id

        experiment_id = ExperimentId(id)
        self._file_management.delete_experiment_folder(experiment_id)