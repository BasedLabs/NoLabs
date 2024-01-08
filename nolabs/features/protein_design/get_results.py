from nolabs.features.protein_design.services.file_management import FileManagement
from nolabs.domain.experiment import ExperimentId
from nolabs.api_models.protein_design import RunProteinDesignResponse


class GetResultsFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    def handle(self, id: str) -> RunProteinDesignResponse:
        assert id

        experiment_id = ExperimentId(id)
        protein_design_data = self._file_management.get_experiment_data(experiment_id)
        return RunProteinDesignResponse(
            pdbContents=protein_design_data,
            errors=[]
        )
