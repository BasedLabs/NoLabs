from typing import Dict, Any

from nolabs.domain.experiment import ExperimentId
from nolabs.api_models.biobuddy import FunctionCall, FunctionParam, FileData, FunctionCallReturnData, RcsbPdbData, \
    RcsbPdbMetaData
from nolabs.modules.biobuddy.functions.base_function import BiobuddyFunction, FunctionParameterDefinition
from nolabs.modules.drug_discovery.services.target_file_management import TargetsFileManagement
from nolabs.infrastructure.settings import Settings

import rcsb_pdb_query_microservice
from rcsb_pdb_query_microservice import DefaultApi as rcsbApiDefaultApi
from rcsb_pdb_query_microservice import ApiClient as rcsbApiClient
from rcsb_pdb_query_microservice import Configuration as rcsbApiConfiguration


class QueryRcsbPdbByNamesFunction(BiobuddyFunction):
    def __init__(self, settings: Settings,
                 targets_file_management: TargetsFileManagement):
        parameters = [
            FunctionParameterDefinition(name="protein_names",
                                        type="array",
                                        required=True,
                                        description="Query RCSB PDB by protein  names. "
                                                    "If a user asks for to pull targets/proteins but does not specify "
                                                    "ids then invoke this"
                                                    "method. Watch out for the plural nouns, i.e. if a user asks to "
                                                    "pull rhodopsins, then query rhodopsin.",
                                        items_type="string"),
            FunctionParameterDefinition(name="max_results",
                                        type="integer",
                                        required=False,
                                        description="Number of proteins to pull. If a user asks to pull exactly all "
                                                    "of the"
                                                    "proteins from a certain query,"
                                                    "then don't add this parameter. If the user asks to pull some "
                                                    "targets/proteins, then set this parameter to 10 by default.")
        ]
        super().__init__("query_rcsb_pdb_by_protein_names", "Query RCSB PDB by protein names.", parameters)
        self._settings = settings
        self._targets_file_management = targets_file_management

    def execute(self, experiment_id: ExperimentId, arguments: Dict[str, Any]) -> FunctionCall:
        protein_names = arguments["protein_names"]

        max_results = None
        if "max_results" in arguments:
            max_results = arguments["max_results"]
        configuration = rcsbApiConfiguration(
            host=self._settings.rcsb_pdb_query_host,
        )
        results = []
        with rcsbApiClient(configuration=configuration) as client:
            api_instance = rcsbApiDefaultApi(client)
            for protein_name in protein_names:
                request = rcsb_pdb_query_microservice.GetFastaFilesBySearchQueryRequest(search_query=protein_name,
                                                                                        max_results=max_results)
                response = api_instance.fetch_fetch_fastas_by_search_query_post(
                    get_fasta_files_by_search_query_request=request)

            data = []
            for result in response.fasta_contents:
                fasta_contents = result.fasta_contents

                data.append(RcsbPdbData(
                    content=fasta_contents,
                    metadata=RcsbPdbMetaData(
                        link=result.link
                    )
                ))

            return FunctionCall(function_name="query_rcsb_pdb_by_protein_names", parameters=[FunctionParam(
                name="protein_names",
                value=protein_names)],
                data=FunctionCallReturnData(files=data))
