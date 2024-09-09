__all__ = [
    'ConformationsJob'
]

from datetime import datetime
from enum import Enum
from typing import List

from mongoengine import ReferenceField, EmbeddedDocument, FloatField, EmbeddedDocumentListField, \
    CASCADE, StringField, IntField, EnumField, BooleanField, BinaryField, DateTimeField

from domain.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.common import Job, Protein, JobInputError


class Integrator(str, Enum):
    langevin = 'LangevinIntegator'
    langevin_middle = 'LangevinMiddleIntegator'
    nose_hoover = 'NoseHooverIntegrator'
    brownian = 'BrownianIntegrator'
    variable_verlet = 'VariableVerletIntegrator'


class ConformationsTimeline(EmbeddedDocument):
    message: str = StringField()
    error: str | None = StringField(required=False)
    created_at: datetime = DateTimeField(default=datetime.utcnow)


class ConformationsJob(Job):
    # region Inputs
    protein: Protein | None = ReferenceField(Protein, reverse_delete_rule=CASCADE, required=True)

    integrator: Integrator = EnumField(Integrator, default=Integrator.langevin)
    total_frames: int = IntField(default=10000)
    temperature_k: float = FloatField(default=309.75)
    take_frame_every: int = IntField(default=1000)
    step_size: float = FloatField(default=0.002)
    replace_non_standard_residues: bool = BooleanField(default=False)
    add_missing_atoms: bool = BooleanField(default=False)
    add_missing_hydrogens: bool = BooleanField(default=False)
    friction_coeff: float = FloatField(default=1.0)
    ignore_missing_atoms: bool = BooleanField(default=False)

    inputs_set: bool = BooleanField(default=False)

    # endregion

    # region Outputs

    md_content: bytes | None = BinaryField(required=False)
    timeline: List[ConformationsTimeline] = EmbeddedDocumentListField(ConformationsTimeline, required=False)

    # endregion

    def result_valid(self) -> bool:
        return not not self.md_content

    def set_protein(self, protein: Protein):
        self.protein = protein

        self.inputs_updated_at = datetime.utcnow()

    def _input_errors(self) -> List[JobInputError]:
        errors = []

        if not self.protein:
            errors.append([
                JobInputError(
                    message='Protein is undefined',
                    error_code=ErrorCodes.invalid_job_input
                )
            ])

        if not self.protein.pdb_content:
            errors.append([
                JobInputError(
                    message='Protein pdb content is empty',
                    error_code=ErrorCodes.invalid_job_input
                )
            ])

        if not self.total_frames or self.total_frames <= 0:
            errors.append(
                JobInputError(
                    message='Total frames are undefined or 0',
                    error_code=ErrorCodes.invalid_job_input
                )
            )

        if not self.temperature_k or self.temperature_k <= 0:
            errors.append(
                JobInputError(
                    message='Temperature is undefined or 0',
                    error_code=ErrorCodes.invalid_job_input
                )
            )

        if not self.take_frame_every or self.take_frame_every <= 0:
            errors.append(
                JobInputError(
                    message='Take frame parameter is undefined or 0',
                    error_code=ErrorCodes.invalid_job_input
                )
            )

        if not self.step_size or self.step_size <= 0:
            errors.append(
                JobInputError(
                    message='Step size is undefined or 0',
                    error_code=ErrorCodes.invalid_job_input
                )
            )

        if not self.friction_coeff or self.friction_coeff < 0:
            errors.append(
                JobInputError(
                    message='Friction coefficient is undefined or 0',
                    error_code=ErrorCodes.invalid_job_input
                )
            )

        if not self.integrator:
            errors.append(
                JobInputError(
                    message='Integrator is undefined',
                    error_code=ErrorCodes.invalid_job_input
                )
            )

        return errors

    def set_input(self,
                  protein: Protein,
                  integrator: Integrator = Integrator.langevin,
                  total_frames: int = 10000,
                  temperature_k: float = 309.75,
                  take_frame_every: int = 1000,
                  step_size: float = 0.002,
                  replace_non_standard_residues: bool = False,
                  add_missing_atoms: bool = False,
                  add_missing_hydrogens: bool = False,
                  friction_coeff: float = 1.0,
                  ignore_missing_atoms: bool = False,
                  ):
        self.timeline = []
        self.md_content = None

        self.protein = protein
        self.total_frames = total_frames
        self.temperature_k = temperature_k
        self.take_frame_every = take_frame_every
        self.step_size = step_size
        self.replace_non_standard_residues = replace_non_standard_residues
        self.add_missing_atoms = add_missing_atoms
        self.add_missing_hydrogens = add_missing_hydrogens
        self.friction_coeff = friction_coeff
        self.ignore_missing_atoms = ignore_missing_atoms
        self.inputs_set = True
        self.integrator = integrator

        self.input_errors(throw=True)

        self.inputs_updated_at = datetime.utcnow()

    def clear_result(self):
        self.md_content = None
        self.timeline = []

    def append_timeline(self, timeline: ConformationsTimeline):
        if not self.protein:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        self.timeline.append(timeline)

    def input_valid(self) -> bool:
        return self.inputs_set

    def set_result(self,
                   protein: Protein,
                   md_content: bytes | str | None):
        if not self.protein:
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Cannot set a result on a job without inputs')

        if not protein:
            raise NoLabsException(ErrorCodes.protein_is_undefined, 'Protein must be provided in job result set')

        if protein != self.protein:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs,
                                  'Protein does not match with protein this job was run')

        if not protein.pdb_content:
            raise NoLabsException(ErrorCodes.protein_pdb_is_empty,
                                  'Cannot set a result on a job with a protein without pdb content')

        if isinstance(md_content, str):
            md_content = md_content.encode('utf-8')

        self.clear_result()
        self.md_content = md_content
