__all__ = ["ExperimentRemovedEventHandler"]

from nolabs.application.use_cases.workflow.data import ExperimentWorkflowRelation
from nolabs.domain.models.common import ExperimentRemovedEvent
from nolabs.seedwork.domain.event_handlers import DomainEventHandler


class ExperimentRemovedEventHandler(DomainEventHandler[ExperimentRemovedEvent]):
    def handle(self, event: ExperimentRemovedEvent):
        relation: ExperimentWorkflowRelation = ExperimentWorkflowRelation.objects(
            experiment=event.experiment
        ).first()

        if relation:
            relation.delete()
