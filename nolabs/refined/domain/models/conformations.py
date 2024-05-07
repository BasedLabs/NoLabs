__all__ = [
    'ConformationsJob'
]

from datetime import datetime
from enum import Enum
from typing import List

from mongoengine import ReferenceField, EmbeddedDocument, FloatField, EmbeddedDocumentListField, \
    CASCADE, StringField, IntField, EnumField, BooleanField, BinaryField, DateTimeField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Job, Protein


class Integrator(str, Enum):
    langevin = 'LangevinIntegator'
    langevin_middle = 'LangevinMiddleIntegator'
    nose_hoover = 'NoseHooverIntegrator'
    brownian = 'BrownianIntegrator'
    variable_verlet = 'VariableVerletIntegrator'


class ConformationsTimeline(EmbeddedDocument):
    message: str = StringField()
    error: str | None = StringField()
    created_at: datetime = DateTimeField(default=datetime.utcnow)

    def __init__(self, message: str, error: str | None, created_at: datetime, *args, **kwargs):
        super().__init__(message=message, error=error, created_at=created_at, *args, **kwargs)


class ConformationsJob(Job):
    protein: Protein | None = ReferenceField(Protein, required=False, reverse_delete_rule=CASCADE)
    md_content: bytes | None = BinaryField(required=False)
    timeline: List[ConformationsTimeline] = EmbeddedDocumentListField(ConformationsTimeline, required=False)

    integrator: Integrator = EnumField(Integrator, default=Integrator.langevin)
    total_frames: int = IntField(default=10000)
    temperature_k: float = FloatField(default=273.15)
    take_frame_every: int = IntField(default=1000)
    step_size: float = FloatField(default=0.002)
    replace_non_standard_residues: bool = BooleanField(default=False)
    add_missing_atoms: bool = BooleanField(default=False)
    add_missing_hydrogens: bool = BooleanField(default=False)
    friction_coeff: float = FloatField(default=1.0)
    ignore_missing_atoms: bool = BooleanField(default=False)

    def set_input(self,
                  protein: Protein,
                  integrator: Integrator = Integrator.langevin,
                  total_frames: int = 10000,
                  temperature_k: float = 273.15,
                  take_frame_every: int = 1000,
                  step_size: float = 0.002,
                  replace_non_standard_residues: bool = False,
                  add_missing_atoms: bool = False,
                  add_missing_hydrogens: bool = False,
                  friction_coeff: float = 1.0,
                  ignore_missing_atoms: bool = False,
                  ):
        if not protein:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if not protein.pdb_content:
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Cannot run conformations on protein without pdb specified')

        if not total_frames or total_frames <= 0:
            raise NoLabsException(ErrorCodes.invalid_job_input, ['Total frames must be greater than 0'])

        if not temperature_k or temperature_k <= 0:
            raise NoLabsException(ErrorCodes.invalid_job_input, ['Temperature K must be greater than 0'])

        if not take_frame_every or take_frame_every <= 0:
            raise NoLabsException(ErrorCodes.invalid_job_input, ['Take frame rate must be greater than 0'])

        if not step_size or step_size <= 0:
            raise NoLabsException(ErrorCodes.invalid_job_input, ['Step size must be greater than 0'])

        if friction_coeff < 0:
            raise NoLabsException(ErrorCodes.invalid_job_input, ['Friction coefficient must be greater than 0'])

        if not integrator:
            raise NoLabsException(ErrorCodes.invalid_job_input, ['You must specify integrator'])

        self.total_frames = total_frames
        self.temperature_k = temperature_k
        self.take_frame_every = take_frame_every
        self.step_size = step_size
        self.replace_non_standard_residues = replace_non_standard_residues
        self.add_missing_atoms = add_missing_atoms
        self.add_missing_hydrogens = add_missing_hydrogens
        self.friction_coeff = friction_coeff
        self.ignore_missing_atoms = ignore_missing_atoms

    def clear_result(self):
        self.md_content = None
        self.timeline = []

    def append_timeline(self, protein: Protein, timeline: ConformationsTimeline):
        if not self.protein:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if protein != self.protein:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs)

        self.timeline.append(timeline)

    def set_result(self,
                   protein: Protein,
                   md_content: bytes | str | None):
        if not self.protein:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if protein != self.protein:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs)

        if not protein.pdb_content:
            raise NoLabsException(ErrorCodes.protein_pdb_is_empty)

        if isinstance(md_content, str):
            md_content = md_content.encode('utf-8')

        self.md_content = md_content
