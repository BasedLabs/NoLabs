from __future__ import annotations

__all__ = [
    'WorkflowDependencies'
]

from typing import Annotated, Type, Dict

from fastapi import Depends

from nolabs.refined.application.use_cases.binding_pockets.workflow import BindingPocketPredictionComponent
from nolabs.refined.application.use_cases.conformations.workflow import ConformationComponent
from nolabs.refined.application.use_cases.diffdock.workflow import DiffDockComponent
from nolabs.refined.application.use_cases.folding.workflow import FoldingComponent
from nolabs.refined.application.use_cases.gene_ontology.workflow import GeneOntologyComponent
from nolabs.refined.application.use_cases.ligands.workflow import LigandsComponent
from nolabs.refined.application.use_cases.localisation.workflow import LocalisationComponent
from nolabs.refined.application.use_cases.msa_generation.workflow import MsaGenerationComponent
from nolabs.refined.application.use_cases.protein_design.workflow import ProteinDesignComponent
from nolabs.refined.application.use_cases.proteins.workflow import ProteinsComponent
from nolabs.refined.application.use_cases.small_molecules_design.workflow import SmallMoleculesDesignLearningComponent
from nolabs.refined.application.use_cases.solubility.workflow import SolubilityComponent
from nolabs.workflow.application.use_cases import CreateWorkflowSchemaFeature, GetWorkflowSchemaFeature, \
    SetWorkflowSchemaFeature, StartWorkflowFeature, DeleteWorkflowSchemaFeature, AllWorkflowSchemasFeature
from nolabs.workflow.component import Component


class WorkflowDependencies:
    @staticmethod
    def available_components() -> Dict[str, Type[Component]]:
        return {
            FoldingComponent.name: FoldingComponent,
            ProteinsComponent.name: ProteinsComponent,
            SolubilityComponent.name: SolubilityComponent,
            SmallMoleculesDesignLearningComponent.name: SmallMoleculesDesignLearningComponent,
            ProteinDesignComponent.name: ProteinDesignComponent,
            MsaGenerationComponent.name: MsaGenerationComponent,
            LocalisationComponent.name: LocalisationComponent,
            LigandsComponent.name: LigandsComponent,
            GeneOntologyComponent.name: GeneOntologyComponent,
            DiffDockComponent.name: DiffDockComponent,
            ConformationComponent.name: ConformationComponent,
            BindingPocketPredictionComponent.name: BindingPocketPredictionComponent
        }

    @staticmethod
    def create_workflow_schema(components: Annotated[
        Dict[str, Type[Component]], Depends(WorkflowDependencies.available_components)]) -> CreateWorkflowSchemaFeature:
        return CreateWorkflowSchemaFeature(
            available_components=components
        )

    @staticmethod
    def delete_workflow_schema() -> DeleteWorkflowSchemaFeature:
        return DeleteWorkflowSchemaFeature()

    @staticmethod
    def get_workflow_schema() -> GetWorkflowSchemaFeature:
        return GetWorkflowSchemaFeature()

    @staticmethod
    def set_workflow_schema(components: Annotated[
        Dict[str, Type[Component]], Depends(WorkflowDependencies.available_components)]) -> SetWorkflowSchemaFeature:
        return SetWorkflowSchemaFeature(
            available_components=components
        )

    @staticmethod
    def start_workflow(components: Annotated[
        Dict[str, Type[Component]], Depends(WorkflowDependencies.available_components)]) -> StartWorkflowFeature:
        return StartWorkflowFeature(
            available_components=components
        )

    @staticmethod
    def all_workflow_schemas():
        return AllWorkflowSchemasFeature()
