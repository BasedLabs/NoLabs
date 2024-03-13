import json

import requests

from chembl_query.api_models import ChEMBLMoleculeRequest, ChEMBLMoleculeResponse, Molecule, DrugIndicationRequest, \
    DrugIndicationResponse

from chembl_query.loggers import logger

__all__ = ['search_chembl_molecules', 'search_drugs_for_condition', 'search_chembl_drugs']


def search_chembl_molecules(request: ChEMBLMoleculeRequest) -> ChEMBLMoleculeResponse:
    base_url = "https://www.ebi.ac.uk/chembl/api/data/molecule"
    # Specify the format=json parameter to ensure JSON is returned
    query_params = {}
    # Start building the query parameters

    filters = request.filters

    if request.limit:
        query_params = {'limit': request.limit, 'format': 'json'}

    # Add filters to the query parameters
    if filters:
        for field, value in filters.items():
            query_params[field] = value

    # Add ordering to the query parameters
    if request.order_by:
        query_params['order_by'] = request.order_by

    search_url = f"{base_url}?{'&'.join(f'{key}={value}' for key, value in query_params.items())}"

    print(search_url)

    response = requests.get(search_url)

    molecules = []
    if response.status_code == 200:
        try:
            data = response.json()
            print(data)
            for entry in data['molecules']:
                chembl_id = entry.get('molecule_chembl_id')
                molecule = fetch_molecule_details(chembl_id)
                if molecule:
                    molecules.append(molecule)

        except requests.exceptions.JSONDecodeError:
            print("Error decoding JSON response.")
    else:
        print(f"Failed to fetch data, status code: {response.status_code}")

    return ChEMBLMoleculeResponse(molecules=molecules)

def search_chembl_drugs(request: ChEMBLMoleculeRequest) -> ChEMBLMoleculeResponse:
    base_url = "https://www.ebi.ac.uk/chembl/api/data/drugs"
    # Specify the format=json parameter to ensure JSON is returned
    query_params = {}
    # Start building the query parameters
    if request.limit:
        query_params = {'limit': request.limit, 'format': 'json'}

    filters = request.filters

    print(filters)

    # Add filters to the query parameters
    if filters:
        for field, value in filters.items():
            query_params[field] = value

    # Add ordering to the query parameters
    if request.order_by:
        query_params['order_by'] = request.order_by

    search_url = f"{base_url}?{'&'.join(f'{key}={value}' for key, value in query_params.items())}"

    response = requests.get(base_url)

    molecules = []
    if response.status_code == 200:
        try:
            data = response.json()
            for entry in data['drugs']:
                chembl_id = entry.get('molecule_chembl_id')
                molecule = fetch_molecule_details(chembl_id)
                if molecule:
                    molecules.append(molecule)

        except requests.exceptions.JSONDecodeError:
            print("Error decoding JSON response.")
    else:
        print(f"Failed to fetch data, status code: {response.status_code}")

    return ChEMBLMoleculeResponse(molecules=molecules)

def search_drugs_for_condition(request: DrugIndicationRequest) -> DrugIndicationResponse:
    base_url = "https://www.ebi.ac.uk/chembl/api/data/drug_indication"

    query_params = {}
    # Start building the query parameters
    if request.limit:
        query_params = {'limit': request.limit, 'format': 'json'}

    # Add condition query if specified (assuming 'condition' is a valid field for filtering)
    if request.condition:
        query_params['condition__icontains'] = request.condition

    # Add filters to the query parameters
    if request.filters:
        for field, value in request.filters.items():
            query_params[field] = value

    # Add ordering to the query parameters
    if request.order_by:
        query_params['order_by'] = request.order_by

    # Build the full search URL with query parameters
    search_url = f"{base_url}?{'&'.join(f'{key}={value}' for key, value in query_params.items())}"

    response = requests.get(search_url)
    drugs = []

    if response.status_code == 200:
        data = response.json()
        for entry in data['drug_indications']:
            chembl_id = entry.get('molecule_chembl_id')
            molecule = fetch_molecule_details(chembl_id)
            if molecule:
                drugs.append(molecule)

        return DrugIndicationResponse(drugs=drugs, total_count=len(drugs), page=1, pages=1)
    else:
        print(f"Failed to fetch data, status code: {response.status_code}")
        return DrugIndicationResponse(drugs=[], total_count=0, page=1, pages=0)


def fetch_molecule_details(chembl_id: str) -> Molecule | None:
    molecule_url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}?format=json"
    response = requests.get(molecule_url)
    if response.status_code == 200:
        data = response.json()
        molecule_type = data.get('molecule_type', 'N/A')
        pref_name = data.get('pref_name')
        synonyms = data.get('molecule_synonyms', [])
        synonym_names = [synonym['molecule_synonym'] for synonym in synonyms]

        smiles = 'N/A'
        if 'molecule_structures' in data and data.get('molecule_structures'):
            smiles = data['molecule_structures'].get('canonical_smiles', 'N/A')

        molecule_link = f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}/"

        if smiles != 'N/A':
            return Molecule(chembl_id=chembl_id, molecule_type=molecule_type, pref_name=pref_name,
                            synonyms=synonym_names, smiles=smiles, link=molecule_link)
    else:
        print(f"Failed to fetch molecule details for {chembl_id}, status code: {response.status_code}")
        return None