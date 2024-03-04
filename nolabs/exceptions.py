from enum import Enum
from typing import List


class ErrorCodes(Enum):
    unknown_exception = 0
    fix_pdb_error = 1
    conformations_update_metadata_error = 2
    amino_acid_solubility_run_error = 3
    no_amino_acids = 4
    experiment_id_not_found = 5
    amino_acid_localisation_run_error = 6
    protein_design_run_error = 7
    protein_design_update_metadata_error = 8
    drug_discovery_folding_error = 9
    folding_method_unknown = 10


class NoLabsException(Exception):
    def __init__(self, messages: List[str], error_code: ErrorCodes):
        self.messages = messages
        self.error_code = error_code
