from enum import Enum
from typing import List, Any


class ErrorCodes(Enum):
    unknown_exception = 0
    conformations_update_metadata_error = 1
    amino_acid_solubility_run_error = 2
    no_amino_acids = 3
    experiment_not_found = 4
    amino_acid_localisation_run_error = 5
    protein_design_run_error = 6
    protein_design_update_metadata_error = 7
    drug_discovery_folding_error = 8
    folding_method_unknown = 9

    job_not_found = 9

    invalid_aa_name = 10
    invalid_aa_content = 11
    invalid_protein_name = 12
    invalid_protein_content = 13
    invalid_protein_id = 14
    invalid_aa_id = 15
    invalid_experiment_id = 16
    invalid_experiment_name = 17
    invalid_job_id = 18

    invalid_root_directory = 19
    invalid_job_input = 20
    invalid_job_name = 21
    invalid_job_type = 22

    invalid_job_context_id = 23
    invalid_job_output = 24
    invalid_job_metadata_id = 25

    invalid_protein_location_probability = 26
    invalid_job_metadata = 27

    invalid_operation = 28

    invalid_job_state = 29

    duplicate_amino_acid = 30

    # ++++ Biobuddy

    biobuddy_error_generating_response = 31
    biobuddy_unexpected_message_type = 32


class NoLabsException(Exception):
    def __init__(self, messages: List[str] | str, error_code: ErrorCodes):
        if isinstance(messages, str):
            self.messages = [messages]
        else:
            self.messages = messages
        self.error_code = error_code

    @staticmethod
    def throw(error_code: ErrorCodes):
        raise NoLabsException('', error_code)


def ensure_not_none(obj: Any):
    if not obj:
        raise ValueError()
