__all__ = []


from nolabs.infrastructure.mongo_connector import mongo_connect
from nolabs.infrastructure.log import logger
from nolabs.application.event_handlers.di import EventHandlersDependencies

_initialized = False

def initialize():
    global _initialized

    if not _initialized:
        logger.info("Initializing app")

        mongo_connect()
        EventHandlersDependencies.inject()

        from nolabs.application.small_molecules_design.workflow import (
            SmallMoleculesDesignLearningComponent
        )
        from nolabs.application.proteins import ProteinsComponent
        from nolabs.application.diffdock import DiffDockComponent
        from nolabs.application.folding import EsmfoldLightComponent
        from nolabs.application.ligands import LigandsComponent

        from nolabs.domain.workflow.component import ComponentTypeFactory

        ComponentTypeFactory.add_type(SmallMoleculesDesignLearningComponent)
        ComponentTypeFactory.add_type(ProteinsComponent)
        ComponentTypeFactory.add_type(LigandsComponent)
        ComponentTypeFactory.add_type(EsmfoldLightComponent)
        ComponentTypeFactory.add_type(DiffDockComponent)

        _initialized = True