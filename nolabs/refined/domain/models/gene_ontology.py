__all__ = [
    'GeneOntologyJob'
]

from typing import List, Dict, Any, Tuple
from uuid import UUID

from mongoengine import ReferenceField, ListField, PULL, DictField, EmbeddedDocument, UUIDField, \
    EmbeddedDocumentListField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Job, Protein


class GeneOntologyJobResult(EmbeddedDocument):
    protein_id: UUID = UUIDField(required=True)
    gene_ontology: Dict[str, Any] = DictField(required=False)


class GeneOntologyJob(Job):
    # region Inputs

    proteins: List[Protein] = ListField(ReferenceField(Protein, required=True, reverse_delete_rule=PULL))

    # endregion
    gene_ontologies: List[GeneOntologyJobResult] = EmbeddedDocumentListField(GeneOntologyJobResult)

    def set_inputs(self, proteins: List[Protein]):
        if not proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        self.gene_ontologies = []
        self.proteins = proteins

    def clear_result(self):
        self.gene_ontologies = []

    def result_valid(self) -> bool:
        return not not self.gene_ontologies

    def set_result(self, result: List[Tuple[Protein, Dict[str, Any]]]):
        if not self.proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Cannot set a result on a job without inputs')

        if not result:
            raise NoLabsException(ErrorCodes.invalid_job_result)

        self.gene_ontologies = []

        for protein, go in result:
            if not [p for p in self.proteins if p.iid == protein.iid]:
                raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs,
                                      'Cannot set result for this protein since it is not found in job inputs')

            self.gene_ontologies.append(
                GeneOntologyJobResult(
                    protein_id=protein.iid.value,
                    gene_ontology=go
                )
            )
