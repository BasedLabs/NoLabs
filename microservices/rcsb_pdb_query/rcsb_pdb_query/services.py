import requests
import json

from rcsb_pdb_query.api_models import GetFastaFilesByIdsRequest, GetFastaFilesResponse, FetchedProtein, \
    GetFastaFilesBySearchQueryRequest

from rcsb_pdb_query.loggers import logger

__all__ = ['fetch_fasta_files_by_ids', 'fetch_protein_entries_by_name']


def fetch_fasta_files_by_ids(request: GetFastaFilesByIdsRequest) -> GetFastaFilesResponse:
    base_url = "https://www.rcsb.org/fasta/entry/"
    pdb_base_link = "https://www.rcsb.org/structure/"
    results = []

    for pdb_id in request.rcsb_pdb_ids:
        fasta_response = requests.get(f"{base_url}{pdb_id}")
        if fasta_response.status_code == 200:
            results.append(FetchedProtein(fasta_response.text, f"{pdb_base_link}{pdb_id}"))
        else:
            print(f"Failed to fetch FASTA for ID: {pdb_id}")

    return GetFastaFilesResponse(results)


def fetch_protein_entries_by_name(request: GetFastaFilesBySearchQueryRequest) -> GetFastaFilesResponse:
    protein_name = request.search_query
    max_results = request.max_results

    request_options = {}

    if max_results is not None:
        request_options = {
            "return_all_hits": False,
            "paginate": {
                "start": 0,  # Assuming you want to start from the first result
                "rows": max_results  # Limit the number of results to max_results
            }
        }
    else:
        request_options = {
            "return_all_hits": True,
        }
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    fasta_base_url = "https://www.rcsb.org/fasta/entry/"
    pdb_base_link = "https://www.rcsb.org/structure/"
    headers = {'Content-Type': 'application/json'}
    query_payload = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.selected_polymer_entity_types",
                        "operator": "exact_match",
                        "value": "Protein (only)"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "struct.title",
                        "operator": "contains_phrase",
                        "value": protein_name
                    }
                }
            ]
        },
        "request_options": request_options,
        "return_type": "entry"
    }

    search_response = requests.post(search_url, headers=headers, data=json.dumps(query_payload))
    fetched_proteins = []

    if search_response.status_code == 200:
        search_results = search_response.json()

        for result in search_results.get("result_set", []):
            pdb_id = result.get("identifier")
            fasta_response = requests.get(f"{fasta_base_url}{pdb_id}")

            if fasta_response.status_code == 200:
                fasta_contents = fasta_response.text
                link = f"{pdb_base_link}{pdb_id}"
                fetched_proteins.append(FetchedProtein(fasta_contents=fasta_contents, link=link))
            else:
                print(f"Failed to fetch FASTA for PDB ID: {pdb_id}")

        return GetFastaFilesResponse(fasta_contents=fetched_proteins)

    return GetFastaFilesResponse(fasta_contents=[])
