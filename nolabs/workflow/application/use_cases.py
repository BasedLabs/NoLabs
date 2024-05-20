import uuid
from collections import Counter
from typing import Optional, List
from uuid import UUID

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Experiment
from nolabs.workflow.workflow_schema import WorkflowSchemaModel, ComponentModel, PropertyModel
from nolabs.workflow.application.mongoengine_models import WorkflowSchemaDbModel
from nolabs.workflow.function import Component
from nolabs.workflow.component import PythonFunction, Function
from nolabs.workflow.properties import Property
from nolabs.workflow.workflow import Workflow


def map_property(property: Property) -> PropertyModel:
    return PropertyModel(
        type=property.type,
        properties={name: map_property(prop) for name, prop in (property.properties if
                                                                property.properties else {}).items()},
        required=property.required,
        description=property.description,
        enum=property.enum,
        const=property.const,
        format=property.format,
        default=property.default,
        example=property.example,
        title=property.title,
        anyOf=[
            map_property(prop)
            for prop in property.anyOf
        ],
        ref=property.ref
    )


class CreateWorkflowSchemaFeature:
    async def handle(self,
                     experiment_id: UUID,
                     components: List[Component]) -> WorkflowSchemaModel:
        experiment: Experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        workflows: WorkflowSchemaDbModel = WorkflowSchemaDbModel.objects(experiment=experiment)

        if workflows:
            return workflows[0]

        components_models: List[ComponentModel] = []

        ids = Counter([comp.id for comp in components])

        for id, count in ids.items():
            if count > 1:
                raise NoLabsException(ErrorCodes.same_component_already_registered, messages=id)

        for component in components:
            function = PythonFunction(
                component=component
            )
            components_models.append(ComponentModel(
                id=component.id,
                title=component.title,
                input={name: map_property(prop) for name, prop in function.input_properties.items()},
                output={name: map_property(prop) for name, prop in function.input_properties.items()}
            ))

        workflow_schema = WorkflowSchemaModel(
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

        return workflow_schema


class GetWorkflowSchemaFeature:
    async def handle(self, experiment_id: UUID) -> Optional[WorkflowSchemaModel]:
        experiment: Experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        db_models: List[WorkflowSchemaDbModel] = WorkflowSchemaDbModel.objects(experiment=experiment)

        if not db_models:
            return None

        return db_models[0].get_workflow_value()


class SetWorkflowSchemaFeature:
    async def handle(self,
                     workflow_schema: WorkflowSchemaModel,
                     components: List[Component]
                     ) -> WorkflowSchemaModel:
        experiment: Experiment = Experiment.objects.with_id(workflow_schema.experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        Workflow.set_schema_errors(workflow_schema=workflow_schema, components=components)

        db_models: List[WorkflowSchemaDbModel] = WorkflowSchemaDbModel.objects(experiment=experiment)

        if not db_models:
            db_model = WorkflowSchemaDbModel.create(
                id=uuid.uuid4(),
                experiment=experiment,
                value=workflow_schema
            )
        else:
            db_model = db_models[0]
            db_model.set_workflow_value(value=workflow_schema)

        db_model.save(cascade=True)

        return workflow_schema


class StartWorkflowFeature:
    async def handle(self, experiment_id: UUID, components: List[Component]):
        experiment: Experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        db_models: List[WorkflowSchemaDbModel] = WorkflowSchemaDbModel.objects(experiment=experiment)

        if not db_models:
            raise NoLabsException(ErrorCodes.workflow_not_found)

        workflow_schema_model = db_models[0].get_workflow_value()

        workflow = Workflow.create_from_schema(
            workflow_schema_model=workflow_schema_model,
            components=components)

        await workflow.execute(terminate=True)


class StopWorkflowFeature:
    async def handle(self, experiment_id: UUID):
        experiment: Experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        db_models: List[WorkflowSchemaDbModel] = WorkflowSchemaDbModel.objects(experiment=experiment)

        if not db_models:
            raise NoLabsException(ErrorCodes.workflow_not_found)

        workflow_schema_model = db_models[0].get_workflow_value()

        workflow = Workflow.create_from_schema(
            workflow_schema_model=workflow_schema_model,
            components=components)

        await workflow.execute(terminate=True)
