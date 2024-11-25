__all__ = ["TargetResponse", "JobResponse","EstimatesResponse", "GetEstimatesRequest","SetupJobRequest"]

import uuid
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
    number_of_designs: int
    dna_length: int
    replicates: int
    report_email: EmailStr
    target_id: str
    swissprot_id: str
    cart_total: int
    session_url: str


class SetupJobRequest(BaseModel):
    job_id: uuid.UUID
    number_of_designs: int = Field(
        ...,
        ge=24,
        le=500,
        description="Every protein counts as one design. The more designs you test the cheaper it becomes."
    )
    dna_length: int = Field(
        ...,
        ge=0,
        le=300,
        description="AA (avg. protein length)"
    )
    replicates: int = Field(
        ...,
        ge=1,
        le=5,
        description="A replicate is a single test for a protein design. The more replicates you select the more confidence you can have in your data."
    )
    report_email: EmailStr = Field(
        ...,
        description="Email address to send the report"
    )
    target_id: str = Field(
        ...,
        description="Target id of the experiment"
    )
    cart_total: int = Field(
        ...,
        ge=0,
        description="Total price of the experiment"
    )
    swissprot_id: str = Field(...)

