from typing import Any, Dict

import external_data_query_microservice

from nolabs.application.use_cases.biobuddy.api_models import (
    FunctionCall, FunctionCallReturnData, FunctionParam, RcsbPdbData,
    RcsbPdbMetaData)
from nolabs.application.use_cases.biobuddy.functions.base_function import (
    BiobuddyFunction, FunctionParameterDefinition)


class QueryRcsbPdbByIdFunction(BiobuddyFunction):
    def __init__(
        self, rcsb_pdb_query_microservice: external_data_query_microservice.DefaultApi
    ):
        parameters = [
            FunctionParameterDefinition(
                name="pdb_ids",
                type="array",
                required=True,
                description="Query RCSB PDB by protein ids. Don't use DNA or RNA ids. If a "
                "user asks for to pull targets/proteins then invoke this method.",
                items_type="string",
            )
        ]
        super().__init__(
            "query_rcsb_pdb_by_id", "Query RCSB PDB by protein ids.", parameters
        )
        self._rcsb_pdb_query_microservice = rcsb_pdb_query_microservice

    def execute(self, arguments: Dict[str, Any]) -> FunctionCall:
        pdb_ids = arguments[self.parameters[0].name]
        request = external_data_query_microservice.GetFastaFilesByIdsRequest(
            rcsb_pdb_ids=pdb_ids
        )
        response = self._rcsb_pdb_query_microservice.fetch_fetch_fastas_by_ids_post(
            get_fasta_files_by_ids_request=request
        )

        data = []
        for result in response.fasta_contents:
            fasta_contents = result.fasta_contents

            data.append(
                RcsbPdbData(
                    content=fasta_contents, metadata=RcsbPdbMetaData(link=result.link)
                )
            )

        return FunctionCall(
            function_name="query_rcsb_pdb_by_id",
            arguments=[FunctionParam(name="pdb_ids", value=pdb_ids)],
            data=FunctionCallReturnData(files=data),
        )
