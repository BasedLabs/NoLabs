from nolabs.seedwork.domain.event_handlers import DomainEventHandler


class ComponentDeletedEventHandler(DomainEventHandler[ExperimentRemovedEvent]):
    def handle(self, event: ExperimentRemovedEvent):
        relation: ExperimentWorkflowRelation = ExperimentWorkflowRelation.objects(
            experiment=event.experiment
        ).first()

        if relation:
            relation.delete()
