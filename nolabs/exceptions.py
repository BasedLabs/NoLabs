from enum import Enum
from typing import List


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

    # ++++ Small molecules design

    job_not_found = 9

    # ---- Small molecules design


class NoLabsException(Exception):
    def __init__(self, messages: List[str], error_code: ErrorCodes):
        self.messages = messages
        self.error_code = error_code
