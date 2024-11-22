from typing import List

from mongoengine import IntField, ListField, ReferenceField, PULL, EmailField, StringField

from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.domain.models.common import Job, Protein, JobInputError


class ProteinAffinityCharacterizationJob(Job):
    number_of_designs: int = IntField(
        min_value=24,
        max_value=500,
        help_message="Every protein counts as one design. The more designs you test the cheaper it becomes.")
    dna_length: int = IntField(min_value=0, max_value=300, help_message="AA (avg. protein length)")
    replicates: int = IntField(
        min_value=1,
        max_value=5,
        help_message="A replicate is a single test for a protein design. The more replicates you select the more confidence you can have in your data.")
    report_email: str = EmailField(required=True, help_message="Email address to send the report")
    proteins: List[Protein] = ListField(
        ReferenceField(Protein, required=True, reverse_delete_rule=PULL)
    )
    target_id: str = StringField(
        help_message="Target id of the experiment"
    )
    cart_total: int = IntField(
        help_message="Total price of the experiment",
        min_value=0,
    )
    session_url: str = StringField(
        help_message="URL of the session where the experiment is located"
    )
    swissprot_id: str = StringField()

    def set_input(self,
                  number_of_designs: int,
                  dna_length: int,
                  replicates: int,
                  report_email: str,
                  target_id: str,
                  cart_total: int,
                  swissprot_id: str
                  ):
        self.number_of_designs = number_of_designs
        self.dna_length = dna_length
        self.replicates = replicates
        self.report_email = report_email
        self.target_id = target_id
        self.cart_total = cart_total
        self.swissprot_id = swissprot_id

    def set_proteins(self, proteins: List[Protein]):
        if not proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input, "Proteins are empty")

        for protein in proteins:
            if not protein.get_fasta():
                raise NoLabsException(ErrorCodes.invalid_job_input, "Protein fasta content is undefined")

        self.proteins = proteins

    def result_valid(self) -> bool:
        return True

    def _input_errors(self) -> List[JobInputError]:
        errors = []

        if not self.number_of_designs:
            errors.append(
                JobInputError(
                    message="You must specify a number of designs. Every protein counts as one design. The more designs you test the cheaper it becomes.",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if not self.proteins or [p for p in self.proteins if not p.get_fasta()]:
            errors.append(
                JobInputError(
                    message="Protein fasta content is undefined or protein array is empty",
                    error_code=ErrorCodes.protein_fasta_is_empty,
                )
            )

        if not self.replicates:
            errors.append(
                JobInputError(
                    message="A replicate is a single test for a protein design. The more replicates you select the more confidence you can have in your data.",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if not self.report_email:
            errors.append(
                JobInputError(
                    message="Email address is not defined",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )
        if not self.target_id:
            errors.append(
                JobInputError(
                    message="Target id is not defined",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )
        if not self.cart_total:
            errors.append(
                JobInputError(
                    message="Total price of the experiment is not defined",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )
        if not self.session_url:
            errors.append(
                JobInputError(
                    message="Session URL is not defined",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )
        if not self.swissprot_id:
            errors.append(
                JobInputError(
                    message="Swissprot id is not defined",
                    error_code=ErrorCodes.invalid_job_input
                )
            )

        return errors
