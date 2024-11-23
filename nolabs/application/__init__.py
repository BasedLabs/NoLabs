__all__ = ["initialize"]

from nolabs.application.blast.workflow import BlastComponent
from nolabs.application.diffdock.workflow import DiffDockComponent
from nolabs.application.proteinmpnn.workflow import ProteinMPNNComponent
from nolabs.application.event_handlers.di import EventHandlersDependencies
from nolabs.application.folding.workflow import EsmfoldLightComponent, EsmfoldComponent
from nolabs.application.ligands.workflow import LigandsComponent
from nolabs.application.proteins.workflow import ProteinsComponent
from nolabs.application.rfdiffusion.workflow import RfDiffusionComponent
from nolabs.workflow.core.component import ComponentTypeFactory


def initialize():
    ComponentTypeFactory.add_type(ProteinsComponent)
    ComponentTypeFactory.add_type(LigandsComponent)
    ComponentTypeFactory.add_type(EsmfoldLightComponent)
    ComponentTypeFactory.add_type(EsmfoldComponent)
    ComponentTypeFactory.add_type(DiffDockComponent)
    ComponentTypeFactory.add_type(ProteinMPNNComponent)
    ComponentTypeFactory.add_type(RfDiffusionComponent)
    ComponentTypeFactory.add_type(BlastComponent)

    EventHandlersDependencies.inject()
