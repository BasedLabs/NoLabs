__all__ = [
    'LocalisationJob'
]

from typing import List
from uuid import UUID

from mongoengine import ReferenceField, ListField, PULL, EmbeddedDocument, FloatField, EmbeddedDocumentListField, \
    UUIDField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Job, Protein, LocalisationProbability


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
    proteins: List[Protein] = ListField(ReferenceField(Protein, required=False, reverse_delete_rule=PULL))
    probabilities: List[LocalisationJobResult] = EmbeddedDocumentListField(LocalisationJobResult)

    def set_proteins(self, proteins: List[Protein]):
        if not proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        self.proteins = proteins

    def clear_result(self):
        self.probabilities = []

    def set_result(self, protein: Protein, localisation: LocalisationProbability):
        if not self.proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if protein not in self.proteins:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs)

        if not protein.fasta_content:
            raise NoLabsException(ErrorCodes.protein_amino_acid_sequence_not_found)

        existing_result = [res for res in self.probabilities if res.protein_id == protein.id]
        if existing_result:
            localisation_result = existing_result[0]
            localisation_result.protein_id = protein.id
            localisation_result.cytosolic = localisation.cytosolic
            localisation_result.mitochondrial = localisation.mitochondrial
            localisation_result.nuclear = localisation.nuclear
            localisation_result.other = localisation.other
            localisation_result.extracellular = localisation.extracellular
        else:
            localisation_result = LocalisationJobResult(
                protein_id=protein.id,
                cytosolic=localisation.cytosolic,
                mitochondrial=localisation.mitochondrial,
                nuclear=localisation.nuclear,
                other=localisation.other,
                extracellular=localisation.extracellular,
            )
            self.probabilities.append(localisation_result)
