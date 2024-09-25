from typing import Annotated

import external_data_query_microservice
from application.biobuddy.functions.query_chembl import QueryChemblFunction
from application.biobuddy.functions.query_chembl_by_disease import \
    QueryChemblByConditionFunction
from application.biobuddy.functions.query_rcsb_pdb import QueryRCSBPDBFunction
from application.biobuddy.functions.query_rcsb_pdb_by_id import \
    QueryRcsbPdbByIdFunction
from fastapi import Depends

from nolabs.infrastructure.di import InfrastructureDependencies


class FunctionDependencies:
    @staticmethod
    def query_chembl(
        chembl_miscroservice: Annotated[
            external_data_query_microservice.DefaultApi,
            Depends(InfrastructureDependencies.external_query_microservice),
        ]
    ) -> QueryChemblFunction:
        return QueryChemblFunction(chembl_microservice=chembl_miscroservice)

    @staticmethod
    def query_chembl_by_disease(
        chembl_miscroservice: Annotated[
            external_data_query_microservice.DefaultApi,
            Depends(InfrastructureDependencies.external_query_microservice),
        ]
    ) -> QueryChemblByConditionFunction:
        return QueryChemblByConditionFunction(chembl_microservice=chembl_miscroservice)

    @staticmethod
    def query_rcsb_pdb_by_id(
        rcsb_pdb_query_microservice: Annotated[
            external_data_query_microservice.DefaultApi,
            Depends(InfrastructureDependencies.external_query_microservice),
        ]
    ) -> QueryRcsbPdbByIdFunction:
        return QueryRcsbPdbByIdFunction(
            rcsb_pdb_query_microservice=rcsb_pdb_query_microservice
        )

    @staticmethod
    def query_rcsb_pdb(
        rcsb_pdb_query_microservice: Annotated[
            external_data_query_microservice.DefaultApi,
            Depends(InfrastructureDependencies.external_query_microservice),
        ]
    ) -> QueryRCSBPDBFunction:
        return QueryRCSBPDBFunction(
            rcsb_pdb_query_microservice=rcsb_pdb_query_microservice
        )
