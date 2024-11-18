from __future__ import annotations

from enum import Enum
from typing import Any, Dict, List, Optional, Union

from external_data_query.mixins import BaseModelMixin
from pydantic import BaseModel
from pydantic.dataclasses import dataclass

# RCSB PDB


@dataclass
class FetchedProtein:
    fasta_contents: str
    link: str


@dataclass
class GetFastaFilesByIdsRequest(BaseModelMixin):
    rcsb_pdb_ids: List[str]
    job_id: Optional[str] = None


@dataclass
class GetFastaFilesResponse(BaseModelMixin):
    fasta_contents: List[FetchedProtein]


@dataclass
class GetFastaFilesBySearchQueryRequest(BaseModelMixin):
    search_query: str
    max_results: int = None
    exact_match: bool = False  # Added to support exact match queries
    job_id: Optional[str] = None


class SequenceType(str, Enum):
    PROTEIN = "protein"
    DNA = "dna"
    RNA = "rna"


class SequenceQueryRequest(BaseModel):
    sequence: str
    sequence_type: SequenceType
    identity_cutoff: Optional[float] = None
    evalue_cutoff: Optional[float] = None
    job_id: Optional[str] = None


@dataclass
class IsJobRunningResponse:
    is_running: bool


class QueryNode(BaseModel):
    pass


class TerminalNode(QueryNode):
    attribute: str
    value: Union[str, int, float, bool, List[str], Dict[str, Any]]
    operator: str
    negation: bool = False


class LogicalNode(QueryNode):
    operator: str
    children: List[QueryNode]


class ComplexQueryRequest(BaseModel):
    query: QueryNode
    max_results: Optional[int] = None
    job_id: Optional[str] = None
