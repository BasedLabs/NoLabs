from enum import Enum
from typing import List, Any

from pydantic.dataclasses import dataclass


@dataclass
class ErrorCode:
    """
    Will be used to represent and describe the error in application
    """
    code: int
    description: str


class ErrorCodes(Enum):
    unknown_exception = ErrorCode(code=0, description="Unknown exception")
    conformations_update_metadata_error = ErrorCode(code=1, description="Error while updating conformations metadata")
    amino_acid_solubility_run_error = ErrorCode(code=2, description="Amino acid solubility run error")
    no_amino_acids = ErrorCode(code=3, description="Amino acids weren't provided")
    experiment_not_found = ErrorCode(code=4, description="Experiment not found")
    amino_acid_localisation_run_error = ErrorCode(code=5, description="Amino acid localisation run error")
    protein_design_run_error = ErrorCode(code=6, description="Protein design run error")
    protein_design_update_metadata_error = ErrorCode(code=7, description="Error while updating protein design metadata")
    drug_discovery_folding_error = ErrorCode(code=8, description="Drug discovery folding error")
    folding_method_unknown = ErrorCode(code=9, description="Folding method is unknown")
    job_not_found = ErrorCode(code=10, description="Job not found")
    invalid_protein_name = ErrorCode(code=13, description="Invalid protein name")
    invalid_protein_content = ErrorCode(code=14, description="Invalid protein content")
    invalid_protein_id = ErrorCode(code=15, description="Invalid protein id")
    invalid_experiment_id = ErrorCode(code=17, description="Invalid experiment id")
    invalid_experiment_name = ErrorCode(code=18, description="Experiment not found")
    invalid_job_id = ErrorCode(code=19, description="Invalid job id")
    invalid_root_directory = ErrorCode(code=20, description="Invalid root directory")
    invalid_job_input = ErrorCode(code=21, description="Invalid job input")
    invalid_job_name = ErrorCode(code=22, description="Invalid job name")
    invalid_job_type = ErrorCode(code=23, description="Invalid job type")
    invalid_job_context_id = ErrorCode(code=24, description="Invalid job context id")
    invalid_job_result = ErrorCode(code=25, description="Invalid job result")
    duplicate_amino_acid = ErrorCode(code=26, description="Duplicate amino acid")
    unknown_localisation_error = ErrorCode(code=27, description="Unknown localisation job error")
    duplicate_protein = ErrorCode(code=28, description="Duplicate protein")
    invalid_localisation_job_output = ErrorCode(code=29, description="Invalid localisation job output")
    invalid_protein_solubility_probability = ErrorCode(code=30, description="Invalid protein solubility probability value")
    no_domain_event_handler = ErrorCode(code=31, description="No domain event handler for this event")
    biobuddy_error_generating_response = ErrorCode(code=32, description="Biobuddy error generating response")
    biobuddy_unexpected_message_type = ErrorCode(code=33, description="Biobuddy unexpected message type")
    invalid_localisation_probability = ErrorCode(code=34, description="Invalid localisation probability value")
    protein_not_found_in_job_inputs = ErrorCode(code=35, description="Protein not found in job inputs")
    protein_amino_acid_sequence_not_found = ErrorCode(code=36, description="Protein amino acid sequence not found")
    invalid_folding_backend = ErrorCode(code=37, description="Invalid folding backend")
    invalid_gene_ontology = ErrorCode(code=38, description="Invalid gene ontology")
    gene_ontology_run_error = ErrorCode(code=39, description="Gene ontology run error")
    invalid_solubility_probability = ErrorCode(code=40, description="Invalid solubility probability")
    unknown_solubility_error = ErrorCode(code=41, description='Unknown error in solubility job')
    unknown_gene_ontology_error = ErrorCode(code=42, description='Unknown error in gene ontology job')
    unknown_folding_error = ErrorCode(code=43, description='Unknown error in folding job')


if len([e.value.code for e in ErrorCodes]) != len(set([e.value.code for e in ErrorCodes])):
    raise ValueError("Invalid ErrorCode initialization")


class NoLabsException(Exception):
    def __init__(self, error_code: ErrorCodes, messages: str | List[str] | None = None):
        self.error_code = error_code.value.code

        if messages:
            if isinstance(messages, str):
                self.messages = [messages]
            else:
                self.messages = messages
        else:
            self.messages = [
                error_code.value.description
            ]
