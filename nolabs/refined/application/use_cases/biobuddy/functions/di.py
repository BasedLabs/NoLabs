from typing import Annotated

from fastapi import Depends

from nolabs.refined.application.use_cases.biobuddy.functions.query_chembl import QueryChemblFunction
from nolabs.refined.application.use_cases.biobuddy.functions.query_chembl_by_disease import \
    QueryChemblByConditionFunction
from nolabs.refined.application.use_cases.biobuddy.functions.query_rcsb_pdb_by_id import QueryRcsbPdbByIdFunction
from nolabs.refined.application.use_cases.biobuddy.functions.query_rcsb_pdb_by_names import QueryRcsbPdbByNamesFunction
from nolabs.refined.infrastructure.di import InfrastructureDependencies

from chembl_query_microservice import DefaultApi as ChemBLDefaultApi


class FunctionDependencies:
    @staticmethod
    def query_chembl(chembl_miscroservice: Annotated[ChemBLDefaultApi, Depends(InfrastructureDependencies.chembl_query_microservice)]
                              ) -> QueryChemblFunction:
        return QueryChemblFunction(chembl_microservice=chembl_miscroservice)

    @staticmethod
    def query_chembl_by_disease(chembl_miscroservice: Annotated[ChemBLDefaultApi, Depends(InfrastructureDependencies.chembl_query_microservice)]
                              ) -> QueryChemblByConditionFunction:
        return QueryChemblByConditionFunction(chembl_microservice=chembl_miscroservice)

    @staticmethod
    def query_rcsb_pdb_by_id(rcsb_pdb_query_microservice: Annotated[ChemBLDefaultApi, Depends(InfrastructureDependencies.rcsb_pdb_query_microservice)]
                              ) -> QueryRcsbPdbByIdFunction:
        return QueryRcsbPdbByIdFunction(rcsb_pdb_query_microservice=rcsb_pdb_query_microservice)

    @staticmethod
    def query_rcsb_pdb_by_name(rcsb_pdb_query_microservice: Annotated[ChemBLDefaultApi, Depends(InfrastructureDependencies.rcsb_pdb_query_microservice)]
                              ) -> QueryRcsbPdbByNamesFunction:
        return QueryRcsbPdbByNamesFunction(rcsb_pdb_query_microservice=rcsb_pdb_query_microservice)