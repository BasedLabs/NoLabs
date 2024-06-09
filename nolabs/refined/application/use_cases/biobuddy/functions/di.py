from typing import Annotated

from fastapi import Depends

from nolabs.refined.application.use_cases.biobuddy.functions.query_chembl import QueryChemblFunction
from nolabs.refined.infrastructure.di import InfrastructureDependencies

from chembl_query_microservice import DefaultApi as ChemBLDefaultApi


class FunctionDependencies:
    @staticmethod
    def query_chembl_function(chembl_miscroservice: Annotated[ChemBLDefaultApi, Depends(InfrastructureDependencies.chembl_query_microservice)]
                              ) -> QueryChemblFunction:
        return QueryChemblFunction(chembl_microservice=chembl_miscroservice)