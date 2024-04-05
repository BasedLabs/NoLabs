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

class QueryRcsbPdbByIdFunction(BiobuddyFunction):
    def __init__(self, settings: Settings,
                 targets_file_management: TargetsFileManagement):
        parameters = [
            FunctionParameterDefinition(name="pdb_ids",
                                        type="array",
                                        required=True,
                                        description="Query RCSB PDB by protein ids. Don't use DNA or RNA ids. If a "
                                                    "user asks for to pull targets/proteins then invoke this method.",
                                        items_type="string")
        ]
        super().__init__("query_rcsb_pdb_by_id", "Query RCSB PDB by protein ids.", parameters)
        self._settings = settings
        self._targets_file_management = targets_file_management

    def execute(self, experiment_id: ExperimentId, arguments: Dict[str, Any]) -> FunctionCall:
        pdb_ids = arguments[self.parameters[0].name]
        print(f"Executing {self.name} with arguments {arguments}")
        configuration = rcsbApiConfiguration(
            host=self._settings.rcsb_pdb_query_host,
        )
        with rcsbApiClient(configuration=configuration) as client:
            api_instance = rcsbApiDefaultApi(client)
            request = rcsb_pdb_query_microservice.GetFastaFilesByIdsRequest(rcsb_pdb_ids=pdb_ids)
            response = api_instance.fetch_fetch_fastas_by_ids_post(get_fasta_files_by_ids_request=request)

            data = []
            for result in response.fasta_contents:
                fasta_contents = result.fasta_contents

                data.append(RcsbPdbData(
                    content=fasta_contents,
                    metadata=RcsbPdbMetaData(
                        link=result.link
                    )
                ))

            return FunctionCall(function_name="query_rcsb_pdb_by_id", parameters=[FunctionParam(name="pdb_ids",
                                                                                                value=pdb_ids)],
                                data=FunctionCallReturnData(files=data))

