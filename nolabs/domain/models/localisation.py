__all__ = [
    'LocalisationJob'
]

import datetime
from typing import List, Tuple
from uuid import UUID

from mongoengine import ReferenceField, ListField, PULL, EmbeddedDocument, FloatField, EmbeddedDocumentListField, \
    UUIDField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.common import Job, Protein, LocalisationProbability, ProteinId, JobInputError


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
        self.proteins = proteins

        self.input_errors(throw=True)

        self.inputs_updated_at = datetime.datetime.utcnow()

    def result_valid(self) -> bool:
        return not not self.probabilities

    def set_result(self, result: List[Tuple[Protein, LocalisationProbability]]):
        self.probabilities = []

        for protein, prob in result:
            if not [p for p in self.proteins if p.iid == protein.iid]:
                raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs,
                                      'This protein is not in input proteins')

            self.probabilities.append(LocalisationJobResult(
                protein_id=protein.iid.value,
                cytosolic=prob.cytosolic,
                mitochondrial=prob.mitochondrial,
                nuclear=prob.nuclear,
                other=prob.other,
                extracellular=prob.extracellular
            ))

        self.input_errors(throw=True)

    def _input_errors(self) -> List[JobInputError]:
        if not self.proteins:
            return [
                JobInputError(
                    message='Proteins are undefined',
                    error_code=ErrorCodes.invalid_job_input
                )
            ]

        for protein in self.proteins:
            if not protein.fasta_content:
                return [JobInputError(
                    message='Protein fasta content is empty',
                    error_code=ErrorCodes.protein_fasta_is_empty
                )]

        return []
