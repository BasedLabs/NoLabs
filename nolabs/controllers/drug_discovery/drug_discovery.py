from typing import Annotated, Dict

from fastapi import APIRouter, Depends



router = APIRouter(
    prefix='/api/v1/drug-discovery',
    tags=['drug-discovery']
)

