__all__ = [
    'MsaGenerationJob'
]

from typing import List
from uuid import UUID

from mongoengine import ReferenceField, ListField, PULL, EmbeddedDocument, FloatField, EmbeddedDocumentListField, \
    UUIDField, CASCADE, StringField, IntField, BinaryField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Job, Protein, LocalisationProbability


class MsaGenerationJob(Job):
    # region Inputs

    protein: Protein | None = ReferenceField(Protein, required=True, reverse_delete_rule=CASCADE)
    msa: bytes | None = BinaryField(required=True)

    # endregion

    contig: str = StringField(required=False)
    number_of_designs: int = IntField(required=False, default=2)
    hotspots: str = StringField(required=False)
    timesteps: int = IntField(required=False, default=50)

    def set_input(self, protein: Protein):
        if not protein:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if not protein.fasta_content:
            raise NoLabsException(ErrorCodes.protein_fasta_is_empty, 'Cannot run msa job on empty fasta')

        self.clear_result()
        self.protein = protein

    def clear_result(self):
        self.msa = None

    def set_result(self,
                   protein: Protein,
                   msa: bytes | str):
        if not protein:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if protein != self.protein:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs)

        if not protein.pdb_content:
            raise NoLabsException(ErrorCodes.protein_pdb_is_empty)

        if not msa:
            raise NoLabsException(ErrorCodes.invalid_msa)

        if isinstance(msa, str):
            msa = msa.encode('utf-8')

        self.msa = msa
