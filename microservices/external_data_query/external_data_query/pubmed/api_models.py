from __future__ import annotations

from pydantic.dataclasses import dataclass
from typing import List, Optional

from external_data_query.mixins import BaseModelMixin

@dataclass
class FetchedArticle:
    title: str
    summary: str
    link: str


@dataclass
class PubMedSearchRequest(BaseModelMixin):
    search_terms: str
    max_results: int
    api_key: Optional[str] = None
    job_id: Optional[str] = None

@dataclass
class PubMedSearchResponse(BaseModelMixin):
    articles: List[FetchedArticle]


@dataclass
class IsJobRunningResponse:
    is_running: bool
