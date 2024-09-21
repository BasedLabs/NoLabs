from enum import Enum
from typing import Any, Dict, List

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
    conformations_update_metadata_error = ErrorCode(
        code=1, description="Error while updating conformations metadata"
    )
    amino_acid_solubility_run_error = ErrorCode(
        code=2, description="Amino acid solubility run error"
    )
    no_amino_acids = ErrorCode(code=3, description="Amino acids weren't provided")
    experiment_not_found = ErrorCode(code=4, description="Experiment not found")
    amino_acid_localisation_run_error = ErrorCode(
        code=5, description="Amino acid localisation run error"
    )
    protein_design_run_error = ErrorCode(code=6, description="Protein design run error")
    protein_design_update_metadata_error = ErrorCode(
        code=7, description="Error while updating protein design metadata"
    )
    drug_discovery_folding_error = ErrorCode(
        code=8, description="Drug discovery folding error"
    )
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
    unknown_localisation_error = ErrorCode(
        code=27, description="Unknown localisation job error"
    )
    duplicate_protein = ErrorCode(code=28, description="Duplicate protein")
    invalid_localisation_job_output = ErrorCode(
        code=29, description="Invalid localisation job output"
    )
    invalid_protein_solubility_probability = ErrorCode(
        code=30, description="Invalid protein solubility probability value"
    )
    no_domain_event_handler = ErrorCode(
        code=31, description="No domain event handler for this event"
    )
    biobuddy_error_generating_response = ErrorCode(
        code=32, description="Biobuddy error generating response"
    )
    biobuddy_unexpected_message_type = ErrorCode(
        code=33, description="Biobuddy unexpected message type"
    )
    invalid_localisation_probability = ErrorCode(
        code=34, description="Invalid localisation probability value"
    )
    protein_not_found_in_job_inputs = ErrorCode(
        code=35, description="Protein not found in job inputs"
    )
    protein_amino_acid_sequence_not_found = ErrorCode(
        code=36, description="Protein amino acid sequence not found"
    )
    invalid_folding_backend = ErrorCode(code=37, description="Invalid folding backend")
    invalid_gene_ontology = ErrorCode(code=38, description="Invalid gene ontology")
    gene_ontology_run_error = ErrorCode(code=39, description="Gene ontology run error")
    invalid_solubility_probability = ErrorCode(
        code=40, description="Invalid solubility probability"
    )
    unknown_solubility_error = ErrorCode(
        code=41, description="Unknown error in solubility job"
    )
    unknown_gene_ontology_error = ErrorCode(
        code=42, description="Unknown error in gene ontology job"
    )
    unknown_folding_error = ErrorCode(
        code=43, description="Unknown error in folding job"
    )
    invalid_ligand_name = ErrorCode(code=44, description="Invalid ligand name")
    invalid_ligand_id = ErrorCode(code=45, description="Invalid ligand id")
    duplicate_ligand = ErrorCode(code=46, description="Duplicate ligand")
    invalid_smiles = ErrorCode(code=47, description="Invalid smiles")
    small_molecules_design_empty_output = ErrorCode(
        code=48, description="Empty output in small molecules design job"
    )
    protein_is_undefined = ErrorCode(code=49, description="Protein is undefined")
    protein_pdb_is_empty = ErrorCode(
        code=50, description="Protein pdb content is empty"
    )
    protein_cannot_be_binder_to_itself = ErrorCode(
        code=51, description="Protein cannot be binder to itself"
    )
    empty_binding_pockets = ErrorCode(code=52, description="Empty binding pockets")
    protein_not_found = ErrorCode(code=53, description="Protein not found")
    invalid_msa = ErrorCode(code=54, description="Invalid msa")
    protein_fasta_is_empty = ErrorCode(code=55, description="Protein fasta is empty")
    binding_pockets_prediction_run_error = ErrorCode(
        code=56, description="Binding pockets prediction run error"
    )
    invalid_drug_likeness_score = ErrorCode(
        code=57, description="Invalid drug likeness score"
    )
    invalid_designed_ligand_score = ErrorCode(
        code=58, description="Invalid designed ligand score"
    )
    invalid_binding_backend = ErrorCode(code=59, description="Invalid binding backend")
    ligand_is_undefined = ErrorCode(code=60, description="Ligand is undefined")
    sdf_content_is_undefined = ErrorCode(
        code=61, description="Sdf content is undefined"
    )
    binding_backend_is_undefined = ErrorCode(
        code=62, description="Binding backend is undefined"
    )
    ligand_not_found = ErrorCode(code=63, description="Ligand not found")
    ligand_initialization_error = ErrorCode(
        code=64, description="Ligand initialization error"
    )
    diffdock_api_error = ErrorCode(
        code=65, description="Diffdock microservice api error"
    )
    protein_msa_is_empty = ErrorCode(code=66, description="Protein msa is empty")
    ligand_smiles_is_empty = ErrorCode(code=67, description="Ligand smiles is empty")
    ligand_not_found_in_job_inputs = ErrorCode(
        code=68, description="Ligand not found in job inputs"
    )
    protein_initialization_error = ErrorCode(
        code=69, description="Protein initialization error"
    )
    fasta_file_is_invalid = ErrorCode(code=70, description="Fasta file is invalid")
    pdb_file_is_invalid = ErrorCode(code=71, description="Pdb file is invalid")
    invalid_ligand_content = ErrorCode(code=72, description="Invalid ligand content")
    smiles_file_is_invalid = ErrorCode(code=73, description="Smiles file is invalid")
    sdf_file_is_invalid = ErrorCode(code=74, description="Sdf file is invalid")
    reinvent_cannot_run_sampling = ErrorCode(code=75, description="Cannot run sampling")
    folding_run_error = ErrorCode(code=76, description="Folding run error")
    workflow_not_found = ErrorCode(code=77, description="Workflow not found")
    same_component_already_registered = ErrorCode(
        code=78, description="Same component was already registered"
    )
    component_has_unmapped_properties = ErrorCode(
        code=79, description="Component has unmapped properties"
    )
    cannot_start_component = ErrorCode(code=80, description="Cannot start component")
    invalid_workflow_schema = ErrorCode(code=81, description="Invalid workflow schema")
    component_not_found = ErrorCode(code=82, description="Component not found")
    component_input_invalid = ErrorCode(
        code=83, description="Component input is invalid"
    )
    flow_run_id_not_found = ErrorCode(code=84, description="Flow run id not found")
    blast_api_error = ErrorCode(code=85, description="Blast api error")

    create_workflow_failed = ErrorCode(code=86, description="Create workflow failed")
    delete_workflow_failed = ErrorCode(code=87, description="Delete workflow failed")
    get_all_workflows_failed = ErrorCode(
        code=88, description="Get all workflows failed"
    )
    get_workflow_schema_failed = ErrorCode(
        code=89, description="Get workflow schema failed"
    )
    update_workflow_schema_failed = ErrorCode(
        code=90, description="Update workflow schema failed"
    )
    start_workflow_schema_failed = ErrorCode(
        code=91, description="Start workflow failed"
    )
    start_component_failed = ErrorCode(code=92, description="Start component failed")
    get_component_state_failed = ErrorCode(
        code=93, description="Get component state failed"
    )
    get_job_state_failed = ErrorCode(code=94, description="Get job state failed")
    job_execution_failed = ErrorCode(code=95, description="Job execution failed")
    get_job_metadata_failed = ErrorCode(code=96, description="Get job metadata failed")
    job_run_not_found = ErrorCode(code=97, description="Job run not found")


if len([e.value.code for e in ErrorCodes]) != len(
    set([e.value.code for e in ErrorCodes])
):
    raise ValueError("Invalid ErrorCode initialization")


class NoLabsException(Exception):
    data: Dict[str, Any] = {}
    message: str

    def __init__(self, error_code: ErrorCodes, message: str | None = None, data: Dict[str, Any] | None = None):
        self.error_code = error_code.value.code

        if message:
            message = message
        else:
            message = error_code.value.description

        self.message = message

        if data is not None:
            self.data = data

        super(Exception, self).__init__(self.message)
