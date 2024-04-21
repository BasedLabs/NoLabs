__all__ = [
    'LocalisationJob',
    'LocalisationProbability',
]

from dataclasses import field
from typing import List

from pydantic.dataclasses import dataclass

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.common.entities import Job, AminoAcidId, JobType, JobState
from nolabs.seedwork.domain.value_objects import ValueObject


@dataclass
class LocalisationProbability(ValueObject):
    amino_acid_id: AminoAcidId
    cytosolic_proteins: float
    mitochondrial_proteins: float
    nuclear_proteins: float
    other_proteins: float
    extracellular_secreted_proteins: float

    def __post_init__(self):
        values = [
            self.cytosolic_proteins,
            self.mitochondrial_proteins,
            self.nuclear_proteins,
            self.other_proteins,
            self.extracellular_secreted_proteins
        ]

        for value in values:
            if not value:
                NoLabsException.throw(ErrorCodes.invalid_protein_location_probability)
            if value < 0 or value > 1.0:
                NoLabsException.throw(ErrorCodes.invalid_protein_location_probability)


@dataclass
class LocalisationJob(Job):
    amino_acid_ids: List[AminoAcidId] = field(default_factory=list)
    output: List[LocalisationProbability] = field(default_factory=list)

    def __post_init__(self):
        super().__post_init__()

        if self.job_type != JobType.LOCALISATION:
            raise NoLabsException('Unsupported job type', ErrorCodes.invalid_job_type)

    def set_input(self, input: List[AminoAcidId]):
        if not input:
            raise NoLabsException('Cannot run job with empty inputs', ErrorCodes.invalid_job_input)

        self.amino_acid_ids = input

        # TODO to add domain event for job running once new UI with async jobs is added
        # instead of set_input and set_state there would be start_job method
        # that will raise this domain event

    def set_state(self, state: JobState):
        if not state:
            raise NoLabsException('Invalid job state', ErrorCodes.invalid_job_state)

        if state == JobState.CREATED and not self.amino_acid_ids:
            raise NoLabsException('Job state cannot be set to started if input is None', ErrorCodes.invalid_job_state)

        self.state = state
