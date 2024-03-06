from typing import Dict, Any

from nolabs.domain.experiment import ExperimentId
from nolabs.features.biobuddy.data_models.message import FunctionCall, FunctionParam
from nolabs.features.biobuddy.functions.base_function import BiobuddyFunction, FunctionParameterDefinition
from nolabs.features.drug_discovery.services.ligand_file_management import LigandsFileManagement
from nolabs.infrastructure.settings import Settings

import chembl_query_microservice
from chembl_query_microservice import DefaultApi as chemblApiDefaultApi
from chembl_query_microservice import ApiClient as chemblApiClient
from chembl_query_microservice import Configuration as chemblApiConfiguration

from nolabs.utils.fasta import create_upload_file_from_string
from nolabs.utils.sdf import smiles_to_sdf_string


class QueryChemblFunction(BiobuddyFunction):
    def __init__(self, settings: Settings,
                 ligands_file_management: LigandsFileManagement):
        parameters = [
            FunctionParameterDefinition(name="search_term",
                                        type="string",
                                        required=True,
                                        description="The search term of the molecule for which we want to find the results")
        ]
        super().__init__("query_chembl", "Query Chembl database by search term to find the molecule",
                         parameters)
        self._settings = settings
        self._ligand_file_management = ligands_file_management

    def execute(self, experiment_id: ExperimentId, arguments: Dict[str, Any]) -> FunctionCall:
        search_term = arguments[self.parameters[0].name]
        configuration = chemblApiConfiguration(
            host=self._settings.chembl_query_host,
        )
        with chemblApiClient(configuration=configuration) as client:
            api_instance = chemblApiDefaultApi(client)
            request = chembl_query_microservice.ChEMBLMoleculeRequest(search_term=search_term)
            molecules = api_instance.query_query_chembl_post(ch_embl_molecule_request=request).molecules

            for molecule in molecules:
                molecule_smiles = molecule.smiles
                sdf_content = smiles_to_sdf_string(molecule_smiles)
                file = create_upload_file_from_string(sdf_content, f"{molecule.chembl_id}.sdf")
                additional_metadata = {
                    "link": molecule.link
                }
                self._ligand_file_management.store_lone_ligand(experiment_id, file, additional_metadata)

            return FunctionCall(function_name='query_chembl', parameters=[FunctionParam(name="search_term",
                                                                                        value=search_term)])

