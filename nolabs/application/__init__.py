__all__ = ["initialize"]

from nolabs.application.diffdock.workflow import DiffDockComponent
from nolabs.application.event_handlers.di import EventHandlersDependencies
from nolabs.application.folding import EsmfoldLightComponent, EsmfoldComponent
from nolabs.application.ligands.workflow import LigandsComponent
from nolabs.application.proteins import ProteinsComponent
from nolabs.workflow.core.component import ComponentTypeFactory


def initialize():
    ComponentTypeFactory.add_type(ProteinsComponent)
    ComponentTypeFactory.add_type(LigandsComponent)
    ComponentTypeFactory.add_type(EsmfoldLightComponent)
    ComponentTypeFactory.add_type(DiffDockComponent)

    EventHandlersDependencies.inject()
