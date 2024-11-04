__all__ = ["initialize"]

from nolabs.application.event_handlers.di import EventHandlersDependencies
from nolabs.application.folding import EsmfoldLightComponent
from nolabs.application.proteins import ProteinsComponent
from nolabs.workflow.core.component import ComponentTypeFactory


def initialize():
    # ComponentTypeFactory.add_type(SmallMoleculesDesignLearningComponent)
    ComponentTypeFactory.add_type(ProteinsComponent)
    # ComponentTypeFactory.add_type(LigandsComponent)
    ComponentTypeFactory.add_type(EsmfoldLightComponent)
    # ComponentTypeFactory.add_type(DiffDockComponent)

    EventHandlersDependencies.inject()
