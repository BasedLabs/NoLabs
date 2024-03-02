import requests

from chembl_query.api_models import ChEMBLMoleculeRequest, ChEMBLMoleculeResponse, Molecule

from chembl_query.loggers import logger

__all__ = ['search_chembl_molecules']


def search_chembl_molecules(request: ChEMBLMoleculeRequest) -> ChEMBLMoleculeResponse:
    base_url = "https://www.ebi.ac.uk/chembl/api/data/molecule"
    # Specify the format=json parameter to ensure JSON is returned
    search_url = f"{base_url}/search?q={request.search_term}&format=json"
    response = requests.get(search_url)

    molecules = []
    if response.status_code == 200:
        try:
            data = response.json()
            for entry in data['molecules']:
                chembl_id = entry.get('molecule_chembl_id')
                molecule_type = entry.get('molecule_type', 'N/A')
                pref_name = entry.get('pref_name')
                synonyms = entry.get('molecule_synonyms', [])
                synonym_names = [synonym['molecule_synonym'] for synonym in synonyms]  # Ensure correct key is used
                print(entry.get('molecule_structures', {}))
                if 'molecule_structures' in entry and entry.get('molecule_structures'):
                    smiles = entry['molecule_structures'].get('canonical_smiles', 'N/A')  # Extract SMILES

                    molecules.append(Molecule(chembl_id=chembl_id, molecule_type=molecule_type, pref_name=pref_name,
                                              synonyms=synonym_names, smiles=smiles))
        except requests.exceptions.JSONDecodeError:
            print("Error decoding JSON response.")
    else:
        print(f"Failed to fetch data, status code: {response.status_code}")

    return ChEMBLMoleculeResponse(molecules=molecules)