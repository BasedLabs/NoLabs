from __future__ import annotations

import dataclasses
from typing import List, Optional, Union, Dict, Any

from pydantic import BaseModel

from external_data_query.mixins import BaseModelMixin

# RCSB PDB

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
    exact_match: bool = False  # Added to support exact match queries
    job_id: Optional[str] = None

@dataclasses.dataclass
class SequenceQueryRequest(BaseModelMixin):
    sequence: str
    sequence_type: str  # "protein", "dna", "rna"
    identity_cutoff: Optional[float] = None
    evalue_cutoff: Optional[float] = None
    job_id: Optional[str] = None

@dataclasses.dataclass
class IsJobRunningResponse:
    is_running: bool

def create_query_node(data: Dict[str, Any]) -> 'QueryNode':
    if 'attribute' in data:
        return TerminalNode(**data)
    elif 'operator' in data and 'children' in data:
        children = [create_query_node(child) for child in data['children']]
        return LogicalNode(operator=data['operator'], children=children)
    else:
        raise ValueError('Invalid input for QueryNode')

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