from typing import Dict, Any

from nolabs.refined.application.use_cases.biobuddy.api_models import FunctionCall, FunctionParam, FunctionCallReturnData, RcsbPdbData, \
    RcsbPdbMetaData
from nolabs.refined.application.use_cases.biobuddy.functions.base_function import BiobuddyFunction, FunctionParameterDefinition
import rcsb_pdb_query_microservice


class QueryRcsbPdbByNamesFunction(BiobuddyFunction):
    def __init__(self, rcsb_pdb_query_microservice: rcsb_pdb_query_microservice.DefaultApi):
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
        self._rcsb_pdb_query_microservice = rcsb_pdb_query_microservice

    def execute(self, arguments: Dict[str, Any]) -> FunctionCall:
        protein_names = arguments["protein_names"]

        max_results = None
        if "max_results" in arguments:
            max_results = arguments["max_results"]

        response = None

        for protein_name in protein_names:
            request = rcsb_pdb_query_microservice.GetFastaFilesBySearchQueryRequest(search_query=protein_name,
                                                                                    max_results=max_results)
            response = self._rcsb_pdb_query_microservice.fetch_fetch_fastas_by_search_query_post(
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

        return FunctionCall(function_name="query_rcsb_pdb_by_protein_names", arguments=[FunctionParam(
            name="protein_names",
            value=protein_names)],
            data=FunctionCallReturnData(files=data))
