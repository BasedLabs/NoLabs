import requests

from rcsb_pdb_query.api_models import GetFastaFilesByIdsRequest, GetFastaFilesByIdsResponse, FetchedProtein

from rcsb_pdb_query.loggers import logger

__all__ = ['fetch_fasta_files_by_ids']

def fetch_fasta_files_by_ids(request: GetFastaFilesByIdsRequest) -> GetFastaFilesByIdsResponse:
    base_url = "https://www.rcsb.org/fasta/entry/"
    pdb_base_link = "https://www.rcsb.org/structure/"
    results = []

    for pdb_id in request.rcsb_pdb_ids:
        fasta_response = requests.get(f"{base_url}{pdb_id}")
        if fasta_response.status_code == 200:
            results.append(FetchedProtein(fasta_response.text, f"{pdb_base_link}{pdb_id}"))
        else:
            print(f"Failed to fetch FASTA for ID: {pdb_id}")

    return GetFastaFilesByIdsResponse(results)