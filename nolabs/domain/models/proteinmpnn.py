# proteinmpnn.py

from typing import List, Optional, Dict
from uuid import uuid4

from mongoengine import (
    CASCADE,
    PULL,
    IntField,
    FloatField,
    BooleanField,
    ListField,
    DictField,
    ReferenceField,
    StringField,
    Document,
    UUIDField,
)

from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.domain.models.common import Job, JobInputError, Protein

class ProteinMPNNResult(Document):
    meta = {'collection': 'protein_mpnn_results'}

    id = UUIDField(primary_key=True, default=uuid4)
    sequence = StringField(required=True)
    fasta_content = StringField(required=True)
    score = FloatField()
    global_score = FloatField()
    T = FloatField()
    sample = IntField()
    seq_recovery = FloatField()
    # Removed the job reference

    def save(self, *args, **kwargs):
        """
        Saves the result to the database.
        """
        super().save(*args, **kwargs)

    def to_dict(self) -> Dict:
        """
        Converts the result to a dictionary for serialization.
        """
        return {
            'id': str(self.id),
            'sequence': self.sequence,
            'fasta_content': self.fasta_content,
            'score': self.score,
            'global_score': self.global_score,
            'T': self.T,
            'sample': self.sample,
            'seq_recovery': self.seq_recovery,
        }

class ProteinMPNNJob(Job):
    # region Inputs
    protein: Protein = ReferenceField(
        Protein, reverse_delete_rule=CASCADE, required=True
    )
    num_seq_per_target: int = IntField(default=2, required=False)
    sampling_temp: float = FloatField(default=0.1, required=False)
    seed: int = IntField(default=37, required=False)
    batch_size: int = IntField(default=1, required=False)
    is_homomer: bool = BooleanField(default=False, required=False)
    chains_to_design: List[str] = ListField(StringField(), required=False)
    fixed_positions: Dict[str, List[int]] = DictField(
        field=ListField(IntField()), required=False
    )
    # endregion

    celery_task_id: Optional[str] = StringField()

    results: List[ProteinMPNNResult] = ListField(
        ReferenceField(ProteinMPNNResult, required=True, reverse_delete_rule=PULL)
    )

    def set_task_id(self, task_id: str):
        """
        Assigns the Celery task ID to the job.
        """
        self.celery_task_id = task_id

    def set_input(
        self,
        protein: Protein,
        num_seq_per_target: int = 2,
        sampling_temp: float = 0.1,
        seed: int = 37,
        batch_size: int = 1,
        is_homomer: bool = False,
        chains_to_design: Optional[List[str]] = None,
        fixed_positions: Optional[Dict[str, List[int]]] = None,
    ):
        """
        Sets the input parameters for the job and validates them.
        """
        self.results = []

        self.protein = protein
        self.num_seq_per_target = num_seq_per_target
        self.sampling_temp = sampling_temp
        self.seed = seed
        self.batch_size = batch_size
        self.is_homomer = is_homomer
        self.chains_to_design = chains_to_design
        self.fixed_positions = fixed_positions

        self.input_errors(throw=True)

    def result_valid(self) -> bool:
        """
        Checks if the job has valid results.
        """
        return bool(self.results)

    def input_valid(self) -> bool:
        """
        Checks if the job's input is valid.
        """
        return bool(self.protein and self.protein.pdb_content)

    def set_result(self, results: List[ProteinMPNNResult]):
        """
        Assigns the results to the job after validation.
        """
        if not results:
            raise NoLabsException(ErrorCodes.invalid_job_result, "Result is empty")

        self.results = results

    def _input_errors(self) -> List[JobInputError]:
        """
        Collects input validation errors.
        """
        errors = []

        if not self.protein:
            errors.append(
                JobInputError(
                    message="Protein is undefined",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if self.protein and not self.protein.pdb_content:
            errors.append(
                JobInputError(
                    message="Protein PDB content is undefined",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if self.num_seq_per_target <= 0:
            errors.append(
                JobInputError(
                    message="Number of sequences per target must be greater than 0",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if self.batch_size <= 0:
            errors.append(
                JobInputError(
                    message="Batch size must be greater than 0",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if self.sampling_temp <= 0:
            errors.append(
                JobInputError(
                    message="Sampling temperature must be greater than 0",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        # Additional validations can be added here

        return errors
