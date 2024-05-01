__all__ = [
    'FoldingJobResult',
    'FoldingJob'
]

from enum import Enum
from typing import List
from uuid import UUID

from mongoengine import ReferenceField, ListField, PULL, EmbeddedDocument, FloatField, EmbeddedDocumentListField, \
    UUIDField, BinaryField, StringField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Job, Protein, LocalisationProbability


class FoldingBackendEnum(str, Enum):
    esmfold = 'esmfold'


class FoldingJobResult(EmbeddedDocument):
    """
    Not a domain object
    Currently used just to keep job outputs
    """
    protein_id: UUID = UUIDField(required=True)
    pdb_content: bytes = BinaryField(required=True)
    backend: str = StringField(required=True)


class FoldingJob(Job):
    proteins: List[Protein] = ListField(ReferenceField(Protein, required=False, reverse_delete_rule=PULL))
    foldings: List[FoldingJobResult] = EmbeddedDocumentListField(FoldingJobResult)

    def set_proteins(self, proteins: List[Protein]):
        if not proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        self.proteins = proteins

    def clear_result(self):
        self.foldings = []

    def set_result(self, protein: Protein, backend: FoldingBackendEnum, pdb_content: bytes):
        if not self.proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if not backend:
            raise NoLabsException(ErrorCodes.invalid_folding_backend)

        if protein not in self.proteins:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs)

        if not protein.fasta_content:
            raise NoLabsException(ErrorCodes.protein_amino_acid_sequence_not_found)

        existing_result = [res for res in self.foldings if res.protein_id == protein.id]
        if existing_result:
            folding_result = existing_result[0]
            folding_result.protein_id = protein.id
            folding_result.pdb_content = pdb_content
            folding_result.backend = backend.value
        else:
            result = FoldingJobResult(
                protein_id=protein.id,
                pdb_content=pdb_content,
                backend=backend.value
            )
            self.foldings.append(result)
