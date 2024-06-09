import json
from typing import Dict, Any

from nolabs.domain.experiment import ExperimentId
from nolabs.api_models.biobuddy import FunctionCall, FunctionParam, ChemBLData, ChemBLMetaData, \
    FunctionCallReturnData
from nolabs.modules.biobuddy.functions.base_function import BiobuddyFunction, FunctionParameterDefinition
from nolabs.modules.drug_discovery.services.ligand_file_management import LigandsFileManagement
from nolabs.infrastructure.settings import Settings

import chembl_query_microservice
from chembl_query_microservice import DefaultApi as chemblApiDefaultApi
from chembl_query_microservice import ApiClient as chemblApiClient
from chembl_query_microservice import Configuration as chemblApiConfiguration

from nolabs.utils.sdf import smiles_to_sdf_string


class QueryChemblByConditionFunction(BiobuddyFunction):
    def __init__(self, settings: Settings,
                 ligands_file_management: LigandsFileManagement):
        parameters = [
            FunctionParameterDefinition(name="filters",
                                        type="string",
                                        required=True,
                                        description="Construct a json filter dictionary to query drugs based on "
                                                    "disease using the ChEMBL API,"
                                                    "construct a filter dictionary by pairing"
                                                    "the API's field names with user-defined values and search "
                                                    "conditions. The fields and their supported conditions"
                                                    "are:\n\n"
                                                    "- `drugind_id` supports: `exact`, `range`, `gt`, `gte`, `lt`, "
                                                    "`lte`, `in`, `isnull`.\n"
                                                    "- `efo_id`, `mesh_heading`, `mesh_id` support: `exact`, "
                                                    "`iexact`, `contains`, `icontains`, `istartswith`,"
                                                    "`startswith`, `endswith`, `iendswith`, `search`, `regex`, "
                                                    "`iregex`, `isnull`, `in`.\n"
                                                    "- `efo_term` (the term of the illness) supports: `exact`, "
                                                    "`iexact`, `contains`, `icontains`, `istartswith`,"
                                                    "`startswith`, `endswith`, `iendswith`, `search`, `regex`, "
                                                    "`iregex`, `isnull`, `in`.\n"
                                                    "- `max_phase_for_ind` supports: `exact`, `range`, `gt`, `gte`, "
                                                    "`lt`, `lte`, `in`, `isnull`.\n\n"
                                                    "Use the format `{field_name__condition: value}` for each filter "
                                                    "entry. For instance, to search for drugs related to a"
                                                    "specific disease with a case-insensitive match, use `{"
                                                    "'efo_term__icontains': 'asthma'}`. Combine multiple fields"
                                                    "and conditions as needed to refine your search based on "
                                                    "disease-related parameters. Remember it's json so remove any "
                                                    "unnecessary single quotations from search terms."),
            FunctionParameterDefinition(name="order",
                                        type="string",
                                        required=True,
                                        description="If a user asks to order the results from querying drugs based on "
                                                    "disease in the ChEMBL API, specify an ordering parameter as a "
                                                    "string"
                                                    "based on the available fields. The ordering parameter determines "
                                                    "how the results are sorted. Available ordering fields"
                                                    "for disease-related queries include:\n\n"
                                                    "- `drugind_id`\n"
                                                    "- `molecule_chembl_id`\n"
                                                    "- `parent_molecule_chembl_id`\n"
                                                    "- `max_phase_for_ind`\n"
                                                    "- `mesh_id`\n"
                                                    "- `mesh_heading`\n"
                                                    "- `efo_id`\n"
                                                    "- `efo_term`\n"
                                                    "- `indication_refs`\n\n"
                                                    "Prefix the field name with a '-' to sort in descending order, "
                                                    "e.g., '-max_phase_for_ind' for sorting by the maximum"
                                                    "phase for indication in descending order. Use the field name as "
                                                    "is for ascending order. Combine multiple fields by"
                                                    "separating them with commas for complex sorting."),
            FunctionParameterDefinition(name="max_molecules",
                                        type="integer",
                                        required=True,
                                        description="Amount of molecules to pull. If a user asks to pull all "
                                                    "molecules set this parameter to 10000. If user does not specify "
                                                    "how many, set it to 20.")
        ]
        super().__init__("query_chembl_by_condition", "Query Chembl database by disease/condition to find the "
                                                      "molecules assocoated"
                                                      "with it",
                         parameters)
        self._settings = settings
        self._ligand_file_management = ligands_file_management

    def execute(self, experiment_id: ExperimentId, arguments: Dict[str, Any]) -> FunctionCall:
        filters = str(arguments[self.parameters[0].name])
        if filters:
            filters = filters.replace("'", '"')
            filters = json.loads(filters)
        order = arguments[self.parameters[1].name]
        max_results = arguments[self.parameters[2].name]
        configuration = chemblApiConfiguration(
            host=self._settings.chembl_query_host,
        )
        with chemblApiClient(configuration=configuration) as client:
            api_instance = chemblApiDefaultApi(client)
            request = None
            if order:
                request = chembl_query_microservice.DrugIndicationRequest(filters=filters,
                                                                          order_by=order,
                                                                          limit=max_results)
            else:
                request = chembl_query_microservice.DrugIndicationRequest(filters=filters,
                                                                          limit=max_results)
            molecules = api_instance.query_query_chembl_by_condition_post(request).drugs

            data = []

            for molecule in molecules:
                molecule_smiles = molecule.smiles
                sdf_content = smiles_to_sdf_string(molecule_smiles)

                molecule.pref_name = molecule.pref_name if molecule.pref_name else molecule.chembl_id
                data.append(ChemBLData(
                    content=sdf_content,
                    metadata=ChemBLMetaData(
                        pref_name=molecule.pref_name,
                        chembl_id=molecule.chembl_id,
                        link=molecule.link
                    )
                ))

            return FunctionCall(function_name='query_chembl_by_condition', parameters=[FunctionParam(name="filters",
                                                                                                     value=filters),
                                                                                       FunctionParam(name="order_by",
                                                                                                     value=order),
                                                                                       FunctionParam(name="max_results",
                                                                                                     value=max_results),
                                                                                       ],
                                data=FunctionCallReturnData(files=data))
