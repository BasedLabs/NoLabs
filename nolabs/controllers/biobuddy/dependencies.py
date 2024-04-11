from typing import Annotated, List

from fastapi import Depends

from nolabs.controllers.common_dependencies import settings_dependency
from nolabs.modules.biobuddy.check_biobuddy_enabled_feature import CheckBioBuddyEnabledFeature
from nolabs.modules.biobuddy.file_management import FileManagement
from nolabs.modules.biobuddy.functions.base_function import BiobuddyFunction
from nolabs.modules.biobuddy.functions.query_chembl import QueryChemblFunction
from nolabs.modules.biobuddy.functions.query_chembl_by_disease import QueryChemblByConditionFunction
from nolabs.modules.biobuddy.functions.query_rcsb_pdb_by_id import QueryRcsbPdbByIdFunction
from nolabs.modules.biobuddy.functions.query_rcsb_pdb_by_names import QueryRcsbPdbByNamesFunction
from nolabs.modules.biobuddy.load_conversation_feature import LoadConversationFeature
from nolabs.modules.biobuddy.send_message_feature import SendMessageFeature
from nolabs.modules.drug_discovery.services.ligand_file_management import LigandsFileManagement
from nolabs.modules.drug_discovery.services.target_file_management import TargetsFileManagement
from nolabs.infrastructure.settings import Settings

def check_biobuddy_enabled_dependency() -> CheckBioBuddyEnabledFeature:
    return CheckBioBuddyEnabledFeature()

def file_management_dependency(settings: Annotated[Settings, Depends(settings_dependency)]) -> FileManagement:
    return FileManagement(settings=settings)


def target_file_management_dependency(
        settings: Annotated[Settings, Depends(settings_dependency)]) -> TargetsFileManagement:
    return TargetsFileManagement(settings=settings)


def ligand_file_management_dependency(settings: Annotated[Settings, Depends(settings_dependency)],
                                      experiments_file_management: Annotated[
                                          FileManagement, Depends(file_management_dependency)],
                                      target_file_management: Annotated[TargetsFileManagement,
                                      Depends(target_file_management_dependency)]) -> LigandsFileManagement:
    return LigandsFileManagement(settings=settings, experiments_file_management=experiments_file_management,
                                 targets_file_management=target_file_management)


def load_conversation_dependency(file_management: Annotated[FileManagement,
Depends(file_management_dependency)]) -> LoadConversationFeature:
    return LoadConversationFeature(file_management=file_management)


def query_rcsb_pdb_by_id_function(settings: Annotated[Settings, Depends(settings_dependency)],
                                             target_file_management: Annotated[TargetsFileManagement,
                                             Depends(target_file_management_dependency)]
                                             ) -> QueryRcsbPdbByIdFunction:
    return QueryRcsbPdbByIdFunction(settings=settings, targets_file_management=target_file_management)


def query_rcsb_pdb_by_names_function(settings: Annotated[Settings, Depends(settings_dependency)],
                                                target_file_management: Annotated[TargetsFileManagement,
                                                Depends(target_file_management_dependency)]
                                                ) -> QueryRcsbPdbByNamesFunction:
    return QueryRcsbPdbByNamesFunction(settings=settings, targets_file_management=target_file_management)


def query_chembl_function(settings: Annotated[Settings, Depends(settings_dependency)],
                                     ligands_file_management: Annotated[LigandsFileManagement,
                                     Depends(ligand_file_management_dependency)]
                                     ) -> QueryChemblFunction:
    return QueryChemblFunction(settings=settings, ligands_file_management=ligands_file_management)

def query_chembl_by_condition_function(settings: Annotated[Settings, Depends(settings_dependency)],
                                     ligands_file_management: Annotated[LigandsFileManagement,
                                     Depends(ligand_file_management_dependency)]
                                     ) -> QueryChemblByConditionFunction:
    return QueryChemblByConditionFunction(settings=settings, ligands_file_management=ligands_file_management)

def send_message_to_drug_discovery_dependency(settings: Annotated[Settings, Depends(settings_dependency)],
                                              file_management: Annotated[FileManagement,
                                              Depends(file_management_dependency)],
                                              target_file_management: Annotated[TargetsFileManagement,
                                              Depends(target_file_management_dependency)],
                                              ligands_file_management: Annotated[LigandsFileManagement,
                                              Depends(ligand_file_management_dependency)]
                                              ) -> SendMessageFeature:
    return SendMessageFeature(settings=settings, file_management=file_management,
                              functions=[query_rcsb_pdb_by_id_function(settings, target_file_management),
                                         query_rcsb_pdb_by_names_function(settings, target_file_management),
                                         query_chembl_function(settings, ligands_file_management),
                                         query_chembl_by_condition_function(settings, ligands_file_management)])
