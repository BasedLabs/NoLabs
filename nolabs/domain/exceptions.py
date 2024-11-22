from enum import Enum
from typing import Any, Dict

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
    experiment_not_found = ErrorCode(code=4, description="Experiment not found")
    protein_design_run_error = ErrorCode(code=6, description="Protein design run error")
    folding_method_unknown = ErrorCode(code=9, description="Folding method is unknown")
    job_not_found = ErrorCode(code=10, description="Job not found")
    invalid_protein_name = ErrorCode(code=13, description="Invalid protein name")
    invalid_protein_content = ErrorCode(code=14, description="Invalid protein content")
    invalid_protein_id = ErrorCode(code=15, description="Invalid protein id")
    invalid_experiment_id = ErrorCode(code=17, description="Invalid experiment id")
    invalid_experiment_name = ErrorCode(code=18, description="Experiment not found")
    invalid_job_id = ErrorCode(code=19, description="Invalid job id")
    invalid_job_input = ErrorCode(code=21, description="Invalid job input")
    invalid_job_name = ErrorCode(code=22, description="Invalid job name")
    invalid_job_result = ErrorCode(code=25, description="Invalid job result")
    no_domain_event_handler = ErrorCode(
        code=31, description="No domain event handler for this event"
    )
    biobuddy_error_generating_response = ErrorCode(
        code=32, description="Biobuddy error generating response"
    )
    invalid_localisation_probability = ErrorCode(
        code=34, description="Invalid localisation probability value"
    )
    protein_not_found_in_job_inputs = ErrorCode(
        code=35, description="Protein not found in job inputs"
    )
    invalid_folding_backend = ErrorCode(code=37, description="Invalid folding backend")
    invalid_gene_ontology = ErrorCode(code=38, description="Invalid gene ontology")
    invalid_solubility_probability = ErrorCode(
        code=40, description="Invalid solubility probability"
    )
    invalid_ligand_name = ErrorCode(code=44, description="Invalid ligand name")
    invalid_ligand_id = ErrorCode(code=45, description="Invalid ligand id")
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
    protein_not_found = ErrorCode(code=53, description="Protein not found")
    invalid_msa = ErrorCode(code=54, description="Invalid msa")
    protein_fasta_is_empty = ErrorCode(code=55, description="Protein fasta is empty")
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
    ligand_not_found = ErrorCode(code=63, description="Ligand not found")
    ligand_initialization_error = ErrorCode(
        code=64, description="Ligand initialization error"
    )
    protein_initialization_error = ErrorCode(
        code=69, description="Protein initialization error"
    )
    fasta_file_is_invalid = ErrorCode(code=70, description="Fasta file is invalid")
    pdb_file_is_invalid = ErrorCode(code=71, description="Pdb file is invalid")
    invalid_ligand_content = ErrorCode(code=72, description="Invalid ligand content")
    smiles_file_is_invalid = ErrorCode(code=73, description="Smiles file is invalid")
    sdf_file_is_invalid = ErrorCode(code=74, description="Sdf file is invalid")
    cannot_start_component = ErrorCode(code=80, description="Cannot start component")
    invalid_workflow_schema = ErrorCode(code=81, description="Invalid workflow schema")
    component_not_found = ErrorCode(code=82, description="Component not found")
    component_input_invalid = ErrorCode(
        code=83, description="Component input is invalid"
    )
    create_workflow_failed = ErrorCode(code=86, description="Create workflow failed")
    delete_workflow_failed = ErrorCode(code=87, description="Delete workflow failed")
    get_workflow_schema_failed = ErrorCode(
        code=89, description="Get workflow schema failed"
    )
    update_workflow_schema_failed = ErrorCode(
        code=90, description="Update workflow schema failed"
    )
    start_workflow_failed = ErrorCode(code=91, description="Start workflow failed")
    start_component_failed = ErrorCode(code=92, description="Start component failed")
    get_component_state_failed = ErrorCode(
        code=93, description="Get component state failed"
    )
    get_job_state_failed = ErrorCode(code=94, description="Get job state failed")
    job_execution_failed = ErrorCode(code=95, description="Job execution failed")
    get_folding_job_failed = ErrorCode(code=98, description="Get folding job failed")
    setup_folding_job_failed = ErrorCode(
        code=99, description="Setup folding job failed"
    )
    run_folding_job_failed = ErrorCode(code=100, description="Run folding job failed")
    workflow_running = ErrorCode(code=103, description="Workflow is running")
    component_running = ErrorCode(code=104, description="Component is running")
    invalid_states_transition = ErrorCode(
        code=110, description="Invalid state transition"
    )
    celery_task_executed_but_task_id_not_found = ErrorCode(
        code=111, description="Celery task executed, but task id not found"
    )
    graph_scheduler_timeout = ErrorCode(code=112, description="Graph scheduler timeout")
    cannot_schedule_node = ErrorCode(code=113, description="Cannot schedule node")
    scheduled_components_are_not_in_graph = ErrorCode(
        code=114, description="Scheduled components are not in graph"
    )
    component_output_invalid = ErrorCode(
        code=115, description="Component output is invalid"
    )
    proteins_are_part_of_another_job = ErrorCode(
        code=116, description="Proteins are used in another job already"
    )
    adaptyv_bio_token_not_set = ErrorCode(
        code=117, description="NOLABS_ADAPTYV_BIO_API_TOKEN not set"
    )
    adaptyv_bio_api_base_not_set = ErrorCode(
        code=118, description="NOLABS_ADAPTYV_BIO_API_BASE not set"
    )
    adaptyv_bio_api_unauthorized = ErrorCode(
        code=119, description="Adaptyv bio api unauthorized exception"
    )
    adaptyv_bio_api_too_many_requests = ErrorCode(
        code=120, description="Adaptyv bio api too many requests"
    )
    adaptyv_bio_api_error = ErrorCode(
        code=121, description="Adaptyv bio api error"
    )
    blast_query_exception = ErrorCode(
        code=122, description="Blast query exception"
    )

if len([e.value.code for e in ErrorCodes]) != len(
    set([e.value.code for e in ErrorCodes])
):
    raise ValueError("Invalid ErrorCode initialization")


class NoLabsException(Exception):
    data: Dict[str, Any] = {}
    message: str

    def __init__(
        self,
        error_code: ErrorCodes,
        message: str | None = None,
        data: Dict[str, Any] | None = None,
    ):
        self.error_code = error_code.value.code

        if message:
            message = message
        else:
            message = error_code.value.description

        self.message = message

        if data is not None:
            self.data = data

        super(Exception, self).__init__(self.message)
