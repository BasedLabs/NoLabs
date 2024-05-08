__all__ = [
    'PocketPredictionJob'
]

from typing import List
from uuid import UUID

from mongoengine import ReferenceField, ListField, PULL, EmbeddedDocument, FloatField, EmbeddedDocumentListField, \
    UUIDField, CASCADE, StringField, IntField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Job, Protein, LocalisationProbability


class PocketPredictionJob(Job):
    # region Inputs

    protein: Protein | None = ReferenceField(Protein, required=True, reverse_delete_rule=CASCADE)

    # endregion

    pocket_ids: List[int] = ListField(IntField())

    def set_input(self,
                  protein: Protein
                  ):
        if not protein:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if not protein.pdb_content:
            raise NoLabsException(ErrorCodes.protein_pdb_is_empty,
                                  'Cannot run pocket prediction job on empty protein pdb')

        self.pocket_ids = []
        self.protein = protein

    def set_result(self,
                   protein: Protein,
                   pocket_ids: List[int]):
        if not self.protein:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if protein != self.protein:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs)

        if not protein.pdb_content:
            raise NoLabsException(ErrorCodes.protein_pdb_is_empty)

        self.pocket_ids = pocket_ids
