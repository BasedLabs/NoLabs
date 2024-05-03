from __future__ import annotations

import dataclasses
from typing import List, Optional

from external_data_query.mixins import BaseModelMixin

@dataclasses.dataclass
class FetchedArticle:
    title: str
    summary: str
    link: str


@dataclasses.dataclass
class PubMedSearchRequest(BaseModelMixin):
    search_terms: str
    max_results: int
    api_key: Optional[str] = None
    job_id: Optional[str] = None

@dataclasses.dataclass
class PubMedSearchResponse(BaseModelMixin):
    articles: List[FetchedArticle]


@dataclasses.dataclass
class IsJobRunningResponse:
    is_running: bool
