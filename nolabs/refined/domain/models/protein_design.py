__all__ = [
    'ProteinDesignJob'
]

from typing import List
from uuid import UUID

from mongoengine import ReferenceField, ListField, PULL, EmbeddedDocument, FloatField, EmbeddedDocumentListField, \
    UUIDField, CASCADE, StringField, IntField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Job, Protein, LocalisationProbability


class ProteinDesignJob(Job):
    # region Inputs

    protein: Protein | None = ReferenceField(Protein, required=True, reverse_delete_rule=CASCADE)
    binders: List[Protein] = ListField(ReferenceField(Protein, required=True, reverse_delete_rule=PULL))

    # endregion

    contig: str = StringField(required=False)
    number_of_designs: int = IntField(required=False, default=2)
    hotspots: str = StringField(required=False)
    timesteps: int = IntField(required=False, default=50)

    def set_input(self,
                  protein: Protein,
                  contig: str,
                  number_of_designs: int,
                  hotspots: str,
                  timesteps: int
                  ):
        if not protein:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if not protein.pdb_content:
            raise NoLabsException(ErrorCodes.protein_pdb_is_empty)

        if not contig:
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Contig cannot be empty')

        if not number_of_designs or number_of_designs <= 0:
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Number of designs must be greater than 0')

        if not timesteps or timesteps <= 0:
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Timesteps must be greater than 0')

        self.binders = []
        self.contig = contig
        self.number_of_designs = number_of_designs
        self.hotspots = hotspots
        self.timesteps = timesteps
        self.protein = protein

    def input_valid(self) -> bool:
        if (not self.protein or not self.contig or not self.protein.pdb_content or not self.number_of_designs or
                self.number_of_designs <= 0 or not self.timesteps or self.timesteps <= 0):
            return False

        return True

    def result_valid(self) -> bool:
        return not not self.binders

    def set_result(self,
                   protein: Protein,
                   binders: List[Protein]):
        self.binders = []

        if not self.protein:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if protein != self.protein:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs)

        if not protein.pdb_content:
            raise NoLabsException(ErrorCodes.protein_pdb_is_empty)

        self.binders = binders
