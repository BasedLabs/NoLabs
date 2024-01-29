from typing import Dict

from nolabs.features.gene_ontology.services.file_management import FileManagement
from nolabs.domain.experiment import ExperimentId
from nolabs.api_models.gene_ontology import RunGeneOntologyResponse, RunGeneOntologyResponseDataNode


class GetResultsFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    def handle(self, id: str) -> RunGeneOntologyResponse:
        assert id

        experiment_id = ExperimentId(id)
        gene_ontology_data = self._file_management.get_result(experiment_id)

        res: Dict[str, RunGeneOntologyResponseDataNode] = {}
        for key in gene_ontology_data.keys():
            res[key] = RunGeneOntologyResponseDataNode(
                name=gene_ontology_data[key].name,
                namespace=gene_ontology_data[key].namespace,
                edges={
                    gene_ontology_data[key].edges
                }
            )

        return RunGeneOntologyResponse(
            data=gene_ontology_data,
            errors=[]
        )
