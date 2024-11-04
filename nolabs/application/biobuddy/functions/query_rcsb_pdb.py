from typing import Any, Dict

from external_data_query_microservice.api.default_api import DefaultApi

from nolabs.application.biobuddy.api_models import (
    FunctionCall,
    FunctionCallReturnData,
    FunctionParam,
    RcsbPdbData,
    RcsbPdbMetaData,
)
from nolabs.application.biobuddy.functions.base_function import (
    BiobuddyFunction,
    FunctionParameterDefinition,
)


class QueryRCSBPDBFunction(BiobuddyFunction):
    def __init__(self, rcsb_pdb_query_microservice: DefaultApi):
        parameters = [
            FunctionParameterDefinition(
                name="query",
                type="string",
                required=True,
                description="A JSON object specifying the conditions for querying the RCSB "
                "PDB database. The JSON structure allows for deeply nested "
                "'and/or' conditions to accommodate complex queries involving "
                "multiple criteria. The JSON must include an 'operator' (e.g., "
                "'and', 'or') at each level, a 'children' array containing "
                "individual or further nested condition objects, "
                "and a 'max_results' key specifying the maximum results to "
                "return. Each condition should include 'attribute', 'value', "
                "'operator', and a 'negation' flag. Example JSON:\n"
                "{\n"
                "  'query': {\n"
                "    'operator': 'and',\n"
                "    'children': [\n"
                "{'attribute': 'struct.title', 'value': 'rhodopsin', 'operator': "
                "'contains_phrase', 'negation': false},\n"
                "{'attribute': 'rcsb_entry_info.selected_polymer_entity_types', "
                "'operator': 'exact_match', 'value': 'Protein (only)', "
                "'negation': false},\n"
                "      {\n"
                "        'operator': 'or',\n"
                "        'children': [\n"
                "{'attribute': 'rcsb_entry_info.assembly_count', 'operator': "
                "'greater', 'value': 1, 'negation': false},\n"
                "{'attribute': 'rcsb_entry_info.molecular_weight', 'operator': "
                "'less', 'value': 5000, 'negation': true}\n"
                "        ]\n"
                "      }\n"
                "    ]\n"
                "  },\n"
                "  'max_results': 10\n"
                "}\n"
                "Supported attribute types and their operators include (please make sure you use only allowed values for the ones where it is specified):\n\n"
                "**Integer Attributes** (supports: equals, greater, less, "
                "greater_or_equal, less_or_equal, range, exists):\n"
                "- 'assembly_count': Number of assemblies defined for the entry.\n"
                "- 'branched_entity_count': Number of distinct branched entities "
                "in the structure entry.\n"
                "- 'cis_peptide_count': Number of cis-peptide linkages per "
                "deposited structure model.\n"
                "- 'deposited_atom_count': Number of heavy atom coordinates "
                "records per deposited structure model.\n"
                "- 'deposited_deuterated_water_count': Number of deuterated water "
                "molecules per deposited structure model.\n"
                "- 'deposited_hydrogen_atom_count': Number of hydrogen atom "
                "coordinates records per deposited structure model.\n"
                "- 'deposited_model_count': Number of model structures deposited.\n"
                "- 'deposited_modeled_polymer_monomer_count': Number of modeled "
                "polymer monomers in the deposited coordinate data.\n"
                "- 'deposited_nonpolymer_entity_instance_count': Number of "
                "non-polymer instances in the deposited data set.\n"
                "- 'deposited_polymer_entity_instance_count': Number of polymer "
                "instances in the deposited data set.\n"
                "- 'deposited_polymer_monomer_count': Number of polymer monomers "
                "in sample entity instances in the deposited data set.\n"
                "- 'deposited_solvent_atom_count': Number of heavy solvent atom "
                "coordinates records per deposited structure model.\n"
                "- 'deposited_unmodeled_polymer_monomer_count': Number of "
                "unmodeled polymer monomers in the deposited coordinate data.\n"
                "- 'disulfide_bond_count': Number of disulfide bonds per "
                "deposited structure model.\n"
                "- 'entity_count': Number of distinct polymer, non-polymer, "
                "branched molecular, and solvent entities per deposited structure "
                "model.\n"
                "- 'experimental_method_count': Number of experimental methods "
                "contributing to the structure determination.\n"
                "- 'inter_mol_covalent_bond_count': Number of intermolecular "
                "covalent bonds.\n"
                "- 'inter_mol_metalic_bond_count': Number of intermolecular "
                "metallic bonds.\n"
                "- 'polymer_entity_count': Number of distinct polymer entities in "
                "the structure entry.\n"
                "- 'polymer_entity_count_DNA': Number of distinct DNA polymer "
                "entities.\n"
                "- 'polymer_entity_count_RNA': Number of distinct RNA polymer "
                "entities.\n"
                "- 'polymer_entity_count_nucleic_acid': Number of distinct "
                "nucleic acid polymer entities (DNA or RNA).\n"
                "- 'polymer_entity_count_nucleic_acid_hybrid': Number of distinct "
                "hybrid nucleic acid polymer entities.\n"
                "- 'polymer_entity_count_protein': Number of distinct protein "
                "polymer entities.\n"
                "- 'polymer_entity_taxonomy_count': Number of distinct taxonomies "
                "represented among the polymer entities in the entry.\n"
                "- 'polymer_monomer_count_maximum': Maximum monomer count of a "
                "polymer entity per deposited structure model.\n"
                "- 'polymer_monomer_count_minimum': Minimum monomer count of a "
                "polymer entity per deposited structure model.\n"
                "- 'solvent_entity_count': Number of distinct solvent entities "
                "per deposited structure model.\n"
                "- 'structure_determination_methodology_priority': Indicates the "
                "priority of the value in the structure determination "
                "methodology.\n\n"
                "**Number Attributes** (supports: equals, greater, less, "
                "greater_or_equal, less_or_equal, range, exists):\n"
                "- 'diffrn_radiation_wavelength_maximum': Maximum radiation "
                "wavelength in angstroms.\n"
                "- 'diffrn_radiation_wavelength_minimum': Minimum radiation "
                "wavelength in angstroms.\n"
                "- 'molecular_weight': Molecular mass (KDa) of polymer and "
                "non-polymer entities in the deposited structure entry.\n"
                "- 'nonpolymer_molecular_weight_maximum': Maximum molecular mass "
                "(KDa) of a non-polymer entity.\n"
                "- 'nonpolymer_molecular_weight_minimum': Minimum molecular mass "
                "(KDa) of a non-polymer entity.\n"
                "- 'polymer_molecular_weight_maximum': Maximum molecular mass ("
                "KDa) of a polymer entity in the deposited structure entry.\n"
                "- 'polymer_molecular_weight_minimum': Minimum molecular mass ("
                "KDa) of a polymer entity in the deposited structure entry.\n"
                "- 'resolution_combined': Combined estimates of experimental "
                "resolution contributing to the refined structural model.\n\n"
                "**String Attributes** (supports: in, exact_match, exists):\n"
                "- 'experimental_method': Category of experimental methods used "
                "to determine the structure.\n"
                "- 'na_polymer_entity_types': Nucleic acid polymer entity type "
                "categories describing the entry.\n"
                "- 'polymer_composition': Categories describing the polymer "
                "entity composition for the entry.\n"
                "- 'selected_polymer_entity_types': DON'T USE IF ASKED TO PULL "
                "DNAs SPECIFICALLY. Selected polymer entity type categories "
                "describing the entry. Allowed values: Nucleic acid (only), "
                "Oligosaccharide (only), Other, Protein (only), Protein/NA, "
                "Protein/Oligosaccharide.\n"
                "- 'software_programs_combined': Combined list of software "
                "programs names reported in connection with the production of "
                "this entry.\n"
                "- 'structure_determination_methodology': Indicates if the "
                "structure was determined using experimental or computational "
                "methods.",
            ),
            FunctionParameterDefinition(
                name="max_results",
                type="integer",
                required=True,
                description="Maximum number of results to retrieve from the query.",
            ),
        ]
        super().__init__(
            "query_rcsb_pdb",
            "Query RCSB PDB database to find entries based on specific conditions",
            parameters,
        )
        self._rcsb_pdb_query_microservice = rcsb_pdb_query_microservice

    def build_query(self, query_conditions, max_results):
        query = {"query": query_conditions, "max_results": max_results}
        return query

    def execute(self, arguments: Dict[str, Any]) -> FunctionCall:
        query_conditions = arguments[self.parameters[0].name]
        max_results = arguments[self.parameters[1].name]
        print(f"Executing {self.name} with arguments {arguments}")

        query = self.build_query(
            query_conditions=query_conditions, max_results=max_results
        )

        response = (
            self._rcsb_pdb_query_microservice.fetch_fetch_fastas_by_complex_query_post(
                body=query
            )
        )

        data = []
        for result in response.fasta_contents:
            print(response.fasta_contents)
            fasta_contents = result.fasta_contents

            data.append(
                RcsbPdbData(
                    content=fasta_contents, metadata=RcsbPdbMetaData(link=result.link)
                )
            )

        return FunctionCall(
            function_name="query_rcsb_pdb",
            arguments=[FunctionParam(name="query_json", value=query)],
            data=FunctionCallReturnData(files=data),
        )
