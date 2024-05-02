import requests
import json

from rcsb_pdb_query.api_models import GetFastaFilesByIdsRequest, GetFastaFilesResponse, FetchedProtein, \
    GetFastaFilesBySearchQueryRequest, AttributeQueryRequest, SequenceQueryRequest, ComplexQueryRequest, QueryNode, \
    TerminalNode, LogicalNode

__all__ = ['fetch_fasta_files_by_ids', 'fetch_protein_entries_by_name',
           'fetch_entries_by_complex_query', 'fetch_proteins_by_sequence']


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
    exact_match = request.exact_match  # Use the new field

    request_options = {}
    if max_results is not None:
        request_options = {
            "return_all_hits": False,
            "paginate": {
                "start": 0,
                "rows": max_results
            }
        }
    else:
        request_options = {"return_all_hits": True}

    operator = "exact_match" if exact_match else "contains_phrase"
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
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
                        "operator": operator,
                        "value": protein_name
                    }
                }
            ]
        },
        "request_options": request_options,
        "return_type": "entry"
    }

    print(query_payload)
    search_response = requests.post(search_url, headers=headers, data=json.dumps(query_payload))
    fetched_proteins = []

    fasta_base_url = "https://www.rcsb.org/fasta/entry/"
    pdb_base_link = "https://www.rcsb.org/structure/"

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


def fetch_entries_by_complex_query(request: ComplexQueryRequest) -> GetFastaFilesResponse:
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    headers = {'Content-Type': 'application/json'}
    fasta_base_url = "https://www.rcsb.org/fasta/entry/"
    pdb_base_link = "https://www.rcsb.org/structure/"

    def build_query(node: QueryNode):
        if isinstance(node, TerminalNode):
            parameters = {
                "attribute": node.attribute,
                "operator": node.operator,
                "value": node.value
            }
            # Only add 'negation' if it is True
            if node.negation:
                parameters["negation"] = node.negation
            return {
                "type": "terminal",
                "service": "text",
                "parameters": parameters
            }
        elif isinstance(node, LogicalNode):
            return {
                "type": "group",
                "logical_operator": node.operator,
                "nodes": [build_query(child) for child in node.children]
            }

    query_payload = {
        "query": build_query(request.query),
        "request_options": {
            "return_all_hits": True if request.max_results is None else False,
            "paginate": {
                "start": 0,
                "rows": request.max_results
            } if request.max_results is not None else {}
        },
        "return_type": "entry"
    }

    print(query_payload)

    response = requests.post(search_url, headers=headers, data=json.dumps(query_payload))
    fetched_proteins = []

    if response.status_code == 200:
        search_results = response.json()

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


def fetch_proteins_by_sequence(request: SequenceQueryRequest) -> GetFastaFilesResponse:
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    headers = {'Content-Type': 'application/json'}
    fasta_base_url = "https://www.rcsb.org/fasta/entry/"
    pdb_base_link = "https://www.rcsb.org/structure/"

    sequence_filters = {}
    if request.identity_cutoff is not None:
        sequence_filters["identity_cutoff"] = request.identity_cutoff
    if request.evalue_cutoff is not None:
        sequence_filters["evalue_cutoff"] = request.evalue_cutoff

    query_payload = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "sequence": request.sequence,
                "sequence_type": request.sequence_type,
                **sequence_filters
            }
        },
        "request_options": {
            "return_all_hits": True
        },
        "return_type": "entry"
    }

    response = requests.post(search_url, headers=headers, data=json.dumps(query_payload))
    fetched_proteins = []

    if response.status_code == 200:
        search_results = response.json()

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