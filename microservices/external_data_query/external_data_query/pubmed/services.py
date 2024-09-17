import os

import certifi
import requests
from Bio import Entrez
from external_data_query.pubmed.api_models import (FetchedArticle,
                                                   PubMedSearchRequest,
                                                   PubMedSearchResponse)

os.environ["SSL_CERT_FILE"] = certifi.where()

__all__ = ["search_pubmed"]


def search_pubmed(request: PubMedSearchRequest) -> PubMedSearchResponse:
    Entrez.email = "example@gmail.com"  # Set the email

    # Perform the search
    handle = Entrez.esearch(
        db="pubmed", term=request.search_terms, retmax=request.max_results
    )
    result = Entrez.read(handle)
    handle.close()

    # Get the list of IDs
    id_list = result["IdList"]

    # Fetch details for each ID
    handle = Entrez.efetch(db="pubmed", id=",".join(id_list), retmode="xml")
    articles = Entrez.read(handle)
    handle.close()

    # Process and store article information
    articles_info = []
    for article in articles["PubmedArticle"]:
        title = article["MedlineCitation"]["Article"]["ArticleTitle"]
        try:
            abstract = article["MedlineCitation"]["Article"]["Abstract"][
                "AbstractText"
            ][0]
        except KeyError:  # Abstract might be missing
            abstract = "No abstract available"
        pmid = article["MedlineCitation"]["PMID"]
        link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

        articles_info.append(FetchedArticle(title=title, summary=abstract, link=link))

    return PubMedSearchResponse(articles=articles_info)
