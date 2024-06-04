from typing import List

from nolabs.workflow.component import Component
from nolabs.workflow.workflow_schema import ComponentModel


class WorkflowTestsMixin:
    def create_workflow(self, components: List[Component]):
        models = []

        for component in components:
            models.append(ComponentModel(
                name=component.name,
                input={name: map_property(prop) for name, prop in component.input_properties.items()},
                output={name: map_property(prop) for name, prop in component.output_properties.items()}
            ))

        workflow_schema = WorkflowSchemaModel(
            workflow_id=uuid.uuid4(),
            experiment_id=experiment_id,
            error=None,
            components=components_models,
            workflow_components=[]
        )

        db_model = WorkflowSchemaDbModel.create(
            id=uuid.uuid4(),
            experiment=experiment,
            value=workflow_schema
        )
        db_model.save(cascade=True)