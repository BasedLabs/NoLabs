from nolabs.api_models.experiment import ExperimentMetadataResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.api_models.gene_ontology import (GetExperimentResponse,
                                             AminoAcidResponse, RunGeneOntologyResponseDataNode)
from nolabs.features.gene_ontology.services.file_management import FileManagement


class GetExperimentFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    def handle(self, id: str) -> GetExperimentResponse:
        assert id

        experiment_id = ExperimentId(id)
        metadata = self._file_management.get_experiment_metadata(experiment_id)
        data = self._file_management.get_experiment_data(experiment_id)
        amino_acids = []
        for amino_acid, obo_nodes in data:
            amino_acids.append(
                AminoAcidResponse(
                    name=amino_acid.name,
                    sequence=amino_acid.sequence,
                    go=[
                        RunGeneOntologyResponseDataNode(name=g.name,
                                                        namespace=g.namespace,
                                                        edges=g.edges) for g in obo_nodes
                    ]
                )
            )
        return GetExperimentResponse(
            amino_acids=amino_acids
        )
