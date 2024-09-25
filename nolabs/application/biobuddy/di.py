__all__ = ["BiobuddyDependencies"]

from typing import Annotated

from biobuddy_microservice import DefaultApi
from fastapi import Depends

from nolabs.application.biobuddy.functions.di import FunctionDependencies
from nolabs.application.biobuddy.functions.query_chembl import QueryChemblFunction
from nolabs.application.biobuddy.functions.query_chembl_by_disease import \
    QueryChemblByConditionFunction
from nolabs.application.biobuddy.functions.query_rcsb_pdb import \
    QueryRCSBPDBFunction
from nolabs.application.biobuddy.functions.query_rcsb_pdb_by_id import \
    QueryRcsbPdbByIdFunction
from nolabs.application.biobuddy.use_cases import (
    CheckBioBuddyEnabledFeature, CreateFunctionCallMessageFeature,
    CreateMessageFeature, EditMessageFeature, GetAvailableFunctionCallsFeature,
    LoadConversationFeature, SendActionQueryFeature)
from nolabs.infrastructure.di import InfrastructureDependencies


class BiobuddyDependencies:
    @staticmethod
    def check_biobuddy_enabled() -> CheckBioBuddyEnabledFeature:
        return CheckBioBuddyEnabledFeature()

    @staticmethod
    def load_conversation() -> LoadConversationFeature:
        return LoadConversationFeature()

    @staticmethod
    def create_message() -> CreateMessageFeature:
        return CreateMessageFeature()

    @staticmethod
    def create_function_call_message() -> CreateFunctionCallMessageFeature:
        return CreateFunctionCallMessageFeature()

    @staticmethod
    def edit_message() -> EditMessageFeature:
        return EditMessageFeature()

    @staticmethod
    def send_query(
        biobuddy: Annotated[
            DefaultApi, Depends(InfrastructureDependencies.biobuddy_microservice)
        ],
        query_chembl: Annotated[
            QueryChemblFunction, Depends(FunctionDependencies.query_chembl)
        ],
        query_chembl_by_disease: Annotated[
            QueryChemblByConditionFunction,
            Depends(FunctionDependencies.query_chembl_by_disease),
        ],
        query_rcsb_pdb_by_id: Annotated[
            QueryRcsbPdbByIdFunction, Depends(FunctionDependencies.query_rcsb_pdb_by_id)
        ],
        query_rcsb_pdb: Annotated[
            QueryRCSBPDBFunction, Depends(FunctionDependencies.query_rcsb_pdb)
        ],
    ) -> SendActionQueryFeature:
        return SendActionQueryFeature(
            biobuddy_microservice=biobuddy,
            functions=[
                query_chembl,
                query_chembl_by_disease,
                query_rcsb_pdb_by_id,
                query_rcsb_pdb,
            ],
        )

    @staticmethod
    def get_tools(
        query_chembl: Annotated[
            QueryChemblFunction, Depends(FunctionDependencies.query_chembl)
        ],
        query_chembl_by_disease: Annotated[
            QueryChemblByConditionFunction,
            Depends(FunctionDependencies.query_chembl_by_disease),
        ],
        query_rcsb_pdb_by_id: Annotated[
            QueryRcsbPdbByIdFunction, Depends(FunctionDependencies.query_rcsb_pdb_by_id)
        ],
        query_rcsb_pdb: Annotated[
            QueryRCSBPDBFunction, Depends(FunctionDependencies.query_rcsb_pdb)
        ],
    ) -> GetAvailableFunctionCallsFeature:
        return GetAvailableFunctionCallsFeature(
            functions=[
                query_chembl,
                query_chembl_by_disease,
                query_rcsb_pdb_by_id,
                query_rcsb_pdb,
            ]
        )
