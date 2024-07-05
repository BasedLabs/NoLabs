import json
from typing import Dict, Any


from nolabs.application.use_cases.biobuddy.api_models import FunctionCall, FunctionParam, ChemBLData, ChemBLMetaData, \
    FunctionCallReturnData
from nolabs.application.use_cases.biobuddy.functions.base_function import BiobuddyFunction, FunctionParameterDefinition

import external_data_query_microservice
from nolabs.utils.sdf import smiles_to_sdf_string


class QueryChemblFunction(BiobuddyFunction):
    def __init__(self, chembl_microservice: external_data_query_microservice.DefaultApi):
        parameters = [
            FunctionParameterDefinition(name="filters",
                                        type="string",
                                        required=True,
                                        description="Construct a json filter dictionary for querying the ChEMBL API by "
                                                    "pairing field names with user-defined values"
                                                    "and search conditions. Here's a guide on which conditions each "
                                                    "field supports:\n\n"
                                                    "- `molecule_type`, `pref_name`, `structure_type`, `usan_stem`, "
                                                    "`usan_stem_definition`, `usan_substem` support:"
                                                    "`iexact`, `icontains`, `istartswith`, `iendswith`, `iregex`, "
                                                    "`isnull`, `in`.\n"
                                                    "- `natural_product`, `oral`, `parenteral`, `polymer_flag`, "
                                                    "`prodrug`, `therapeutic_flag`, `topical`,"
                                                    "`withdrawn_flag` support: `exact`, `isnull`.\n"
                                                    "- `usan_year` supports: `exact`, `range`, `gt`, `gte`, `lt`, "
                                                    "`lte`, `in`, `isnull`.\n\n"
                                                    "Use the format `{field_name__condition: value}` for each entry. "
                                                    "For a simple search for 'aspirin' by its"
                                                    "preferred name, use `{'pref_name__iexact': 'aspirin'}`. Combine "
                                                    "multiple fields and conditions as needed"
                                                    "to refine your search. max_phase signifies the phase of the drug "
                                                    "(i.e. approved/pre-clinical etc.)"),
            FunctionParameterDefinition(name="order",
                                        type="string",
                                        required=True,
                                        description="If a user asks to order the results from the ChEMBL API, specify "
                                                    "an ordering"
                                                    "parameter as a string based on the available fields. Otherwise "
                                                    "return an empty string. \n"
                                                    "The ordering parameter determines how the results are sorted. "
                                                    "Available ordering fields include:\n\n"
                                                    "- `availability_type`\n"
                                                    "- `biotherapeutic`\n"
                                                    "- `black_box_warning`\n"
                                                    "- `chebi_par_id`\n"
                                                    "- `chirality`\n"
                                                    "- `dosed_ingredient`\n"
                                                    "- `first_approval`\n"
                                                    "- `first_in_class`\n"
                                                    "- `indication_class`\n"
                                                    "- `inorganic_flag`\n"
                                                    "- `helm_notation`\n"
                                                    "- `max_phase`\n"
                                                    "- `molecule_chembl_id`\n"
                                                    "- `molecule_hierarchy`\n"
                                                    "- `molecule_properties`\n"
                                                    "- `molecule_structures`\n"
                                                    "- `molecule_type`\n"
                                                    "- `natural_product`\n"
                                                    "- `oral`\n"
                                                    "- `parenteral`\n"
                                                    "- `polymer_flag`\n"
                                                    "- `pref_name`\n"
                                                    "- `prodrug`\n"
                                                    "- `structure_type`\n"
                                                    "- `therapeutic_flag`\n"
                                                    "- `topical`\n"
                                                    "- `usan_stem`\n"
                                                    "- `usan_stem_definition`\n"
                                                    "- `usan_substem`\n"
                                                    "- `usan_year`\n"
                                                    "- `withdrawn_flag`\n"
                                                    "- `chemical_probe`\n\n"
                                                    "Prefix the field name with a '-' to sort in descending order, "
                                                    "e.g., '-first_approval' for sorting by first approval date in "
                                                    "descending order."
                                                    "Use the field name as is for ascending order. Combine multiple "
                                                    "fields by separating them with commas for complex sorting."),
            FunctionParameterDefinition(name="max_molecules",
                                        type="integer",
                                        required=True,
                                        description="Amount of molecules to pull. If a user asks to pull all "
                                                    "molecules set this parameter to 10000. If user does not specify "
                                                    "how many, set it to 20.")
        ]
        self._chembl_microservice = chembl_microservice
        super().__init__("query_chembl", "Query Chembl database by search term to find the molecule",
                         parameters)

    def execute(self, arguments: Dict[str, Any]) -> FunctionCall:
        filters = str(arguments[self.parameters[0].name])
        if filters:
            filters = filters.replace("'", '"')
            filters = json.loads(filters)
        else:
            filters = None
        order = arguments[self.parameters[1].name]
        max_results = arguments[self.parameters[2].name]
        request = None
        if order:
            request = external_data_query_microservice.ChEMBLMoleculeRequest(filters=filters,
                                                                      order_by=order,
                                                                      limit=max_results)
        else:
            request = external_data_query_microservice.ChEMBLMoleculeRequest(filters=filters,
                                                                      limit=max_results)

        print(request)
        molecules = self._chembl_microservice.query_query_chembl_post(ch_embl_molecule_request=request).molecules

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

        return FunctionCall(function_name='query_chembl', arguments=[FunctionParam(name="filters",
                                                                                    value=filters),
                                                                      FunctionParam(name="order_by",
                                                                                    value=order),
                                                                      FunctionParam(name="max_results",
                                                                                    value=max_results),
                                                                      ],
                            data=FunctionCallReturnData(files=data))
