__all__ = ["TargetResponse", "JobResponse","EstimatesResponse", "GetEstimatesRequest","SetupJobRequest"]

import uuid
from typing import Optional
from uuid import UUID

from pydantic import BaseModel, Field, EmailStr


class GetTargetsRequest(BaseModel):
    search_query: str


class TargetResponse(BaseModel):
    id: str
    name: str
    description: str
    swissprot_id: str


class EstimatesResponse(BaseModel):
    total_price: int
    turnaround_time: int


class GetEstimatesRequest(BaseModel):
    job_id: uuid.UUID


class JobResponse(BaseModel):
    job_id: UUID
    job_name: str
    number_of_designs: Optional[int]
    dna_length: Optional[int]
    replicates: Optional[int]
    report_email: Optional[EmailStr]
    target_id: Optional[str]
    swissprot_id: Optional[str]
    cart_total: Optional[int]
    session_url: Optional[str]
    submitted: bool = False


class SetupJobRequest(BaseModel):
    job_id: uuid.UUID
    number_of_designs: Optional[int] = Field(
        None,
        ge=24,
        le=500,
        description="Every protein counts as one design. The more designs you test the cheaper it becomes."
    )
    dna_length: Optional[int] = Field(
        None,
        ge=0,
        le=300,
        description="AA (avg. protein length)"
    )
    replicates: Optional[int] = Field(
        None,
        ge=1,
        le=5,
        description="A replicate is a single test for a protein design. The more replicates you select the more confidence you can have in your data."
    )
    report_email: Optional[EmailStr] = Field(
        None,
        description="Email address to send the report"
    )
    target_id: Optional[str] = Field(
        None,
        description="Target id of the experiment"
    )
    cart_total: Optional[int] = Field(
        None,
        ge=0,
        description="Total price of the experiment"
    )
    swissprot_id: Optional[str] = Field(None)

