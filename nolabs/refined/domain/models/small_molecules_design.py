__all__ = [
    'SmallMoleculesDesignJob'
]

from enum import Enum
from typing import List
from uuid import UUID

from mongoengine import ReferenceField, ListField, PULL, EmbeddedDocument, FloatField, EmbeddedDocumentListField, \
    UUIDField, BinaryField, StringField, CASCADE, IntField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Job, Protein, LocalisationProbability, SolubleProbability, Ligand


class SmallMoleculesDesignJob(Job):
    protein: Protein = ReferenceField(Protein, required=False, reverse_delete_rule=CASCADE)

    center_x: float = FloatField(default=0.0)
    center_y: float = FloatField(default=0.0)
    center_z: float = FloatField(default=0.0)
    size_x: float = FloatField(default=5.0)
    size_y: float = FloatField(default=5.0)
    size_z: float = FloatField(default=5.0)
    batch_size: int = IntField(default=128),
    minscore: float = FloatField(default=0.4),
    epochs: int = IntField(default=50)

    sampling_size: int = IntField(default=5)

    ligands: List[Ligand] = ListField(ReferenceField(Ligand, required=False, reverse_delete_rule=PULL))

    def change_sampling_size(self, sampling_size: int):
        if not sampling_size or sampling_size <= 0:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        self.sampling_size = sampling_size

    def set_inputs(self,
                   protein: Protein,
                   center_x: float,
                   center_y: float,
                   center_z: float,
                   size_x: float,
                   size_y: float,
                   size_z: float,
                   batch_size: int,
                   minscore: float,
                   epochs: int):
        if not protein:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if not protein.pdb_content:
            raise NoLabsException(ErrorCodes.protein_pdb_is_empty)

        self.protein = protein

        if center_x is None:
            raise NoLabsException(ErrorCodes.invalid_job_input, ['Center x is required'])

        if center_y is None:
            raise NoLabsException(ErrorCodes.invalid_job_input, ['Center y is required'])

        if center_z is None:
            raise NoLabsException(ErrorCodes.invalid_job_input, ['Center z is required'])

        if size_x is None:
            raise NoLabsException(ErrorCodes.invalid_job_input, ['Size x is required'])

        if size_y is None:
            raise NoLabsException(ErrorCodes.invalid_job_input, ['Size y is required'])

        if size_z is None:
            raise NoLabsException(ErrorCodes.invalid_job_input, ['Size z is required'])

        if batch_size is None or batch_size <= 0:
            raise NoLabsException(ErrorCodes.invalid_job_input, ['Batch size cannot be empty or negative number or 0'])

        if minscore is None or minscore < 0:
            raise NoLabsException(ErrorCodes.invalid_job_input,
                                  ['Minscore cannot be empty or be a negative number or 0'])

        if epochs is None or epochs < 0:
            raise NoLabsException(ErrorCodes.invalid_job_input,
                                  ['Number of epochs cannot be empty or be a negative number'])

        if self.center_x != center_x or \
                self.center_y != center_y or \
                self.center_z != center_z or \
                self.size_x != size_x or \
                self.size_y != size_y or \
                self.size_z != size_z or \
                self.batch_size != batch_size or \
                self.minscore != minscore or \
                self.epochs != epochs:
            self.ligands = []
            self.protein = protein
            self.center_x = center_x
            self.center_y = center_y
            self.center_z = center_z
            self.size_x = size_x
            self.size_y = size_y
            self.size_z = size_z
            self.batch_size = batch_size
            self.minscore = minscore
            self.epochs = epochs

    def set_result(self, protein: Protein, ligands: List[Ligand]):
        if not ligands:
            raise NoLabsException(ErrorCodes.small_molecules_design_empty_output)

        if not self.protein:
            raise NoLabsException(ErrorCodes.protein_is_undefined)

        if protein != protein:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs)

        self.ligands = ligands
