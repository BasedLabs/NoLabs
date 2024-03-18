from __future__ import annotations

import dataclasses
from typing import List, Optional

from rcsb_pdb_query.mixins import BaseModelMixin

@dataclasses.dataclass
class FetchedProtein:
    fasta_contents: str
    link: str

@dataclasses.dataclass
class GetFastaFilesByIdsRequest(BaseModelMixin):
    rcsb_pdb_ids: List[str]
    job_id: Optional[str] = None

@dataclasses.dataclass
class GetFastaFilesResponse(BaseModelMixin):
    fasta_contents: List[FetchedProtein]

@dataclasses.dataclass
class GetFastaFilesBySearchQueryRequest(BaseModelMixin):
    search_query: str
    max_results: int = None
    job_id: Optional[str] = None


@dataclasses.dataclass
class IsJobRunningResponse:
    is_running: bool
