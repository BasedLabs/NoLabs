__all__ = []

from nolabs.application.event_handlers.di import EventHandlersDependencies

from nolabs.application.proteins import ProteinsComponent
# from nolabs.application.diffdock import DiffDockComponent
from nolabs.application.folding import EsmfoldLightComponent
# from nolabs.application.ligands import LigandsComponent

from nolabs.workflow.component import ComponentTypeFactory

# ComponentTypeFactory.add_type(SmallMoleculesDesignLearningComponent)
ComponentTypeFactory.add_type(ProteinsComponent)
# ComponentTypeFactory.add_type(LigandsComponent)
ComponentTypeFactory.add_type(EsmfoldLightComponent)
# ComponentTypeFactory.add_type(DiffDockComponent)

EventHandlersDependencies.inject()
