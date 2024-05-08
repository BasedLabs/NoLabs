__all__ = [
    'LocalisationJob'
]

from typing import List, Tuple
from uuid import UUID

from mongoengine import ReferenceField, ListField, PULL, EmbeddedDocument, FloatField, EmbeddedDocumentListField, \
    UUIDField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Job, Protein, LocalisationProbability, ProteinId


class LocalisationJobResult(EmbeddedDocument):
    """
    Not a domain object
    Currently used just to keep job outputs
    """
    protein_id: UUID = UUIDField(required=True)
    cytosolic: float = FloatField(required=True)
    mitochondrial: float = FloatField(required=True)
    nuclear: float = FloatField(required=True)
    other: float = FloatField(required=True)
    extracellular: float = FloatField(required=True)


class LocalisationJob(Job):
    # region Inputs

    proteins: List[Protein] = ListField(ReferenceField(Protein, required=True, reverse_delete_rule=PULL))

    # endregion

    probabilities: List[LocalisationJobResult] = EmbeddedDocumentListField(LocalisationJobResult)

    def set_input(self, proteins: List[Protein]):
        if not proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        for protein in proteins:
            if not protein.get_fasta():
                raise NoLabsException(ErrorCodes.protein_fasta_is_empty)

        self.proteins = proteins

    def set_result(self, result: List[Tuple[Protein, LocalisationProbability]]):
        if not self.proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if not result:
            raise NoLabsException(ErrorCodes.invalid_job_result)

        self.probabilities = []

        for protein, prob in result:
            if not [p for p in self.proteins if p.iid == protein.iid]:
                raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs, 'This protein is not in input proteins')

            self.probabilities.append(LocalisationJobResult(
                protein_id=protein.iid.value,
                cytosolic=prob.cytosolic,
                mitochondrial=prob.mitochondrial,
                nuclear=prob.nuclear,
                other=prob.other,
                extracellular=prob.extracellular
            ))
