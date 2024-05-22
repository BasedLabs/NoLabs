import uuid
from collections import Counter
from typing import Optional, List
from uuid import UUID

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Experiment
from nolabs.workflow.workflow_schema import WorkflowSchemaModel, ComponentModel, PropertyModel
from nolabs.workflow.application.mongoengine_models import WorkflowSchemaDbModel
from nolabs.workflow.function import PythonFunction
from nolabs.workflow.component import PythonComponent
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
                     functions: List[PythonFunction]) -> WorkflowSchemaModel:
        experiment: Experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        workflows: WorkflowSchemaDbModel = WorkflowSchemaDbModel.objects(experiment=experiment)

        if workflows:
            return workflows[0]

        components_models: List[ComponentModel] = []

        ids = Counter([func.id for func in functions])

        for id, count in ids.items():
            if count > 1:
                raise NoLabsException(ErrorCodes.same_component_already_registered, messages=str(id))

        for function in functions:
            component = PythonComponent(
                function=function
            )
            components_models.append(ComponentModel(
                id=component.id,
                name=component.name,
                input={name: map_property(prop) for name, prop in component.input_properties.items()},
                output={name: map_property(prop) for name, prop in component.output_properties.items()}
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
                     functions: List[PythonFunction]
                     ) -> WorkflowSchemaModel:
        experiment: Experiment = Experiment.objects.with_id(workflow_schema.experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        Workflow.set_schema_errors(workflow_schema=workflow_schema, functions=functions)

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
    async def handle(self, experiment_id: UUID, functions: List[PythonFunction]):
        experiment: Experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        db_models: List[WorkflowSchemaDbModel] = WorkflowSchemaDbModel.objects(experiment=experiment)

        if not db_models:
            raise NoLabsException(ErrorCodes.workflow_not_found)

        workflow_schema_model = db_models[0].get_workflow_value()

        workflow = Workflow.create_from_schema(
            workflow_schema_model=workflow_schema_model,
            functions=functions)

        await workflow.execute(terminate=True)


class StopWorkflowFeature:
    async def handle(self, experiment_id: UUID, functions: List[PythonFunction]):
        experiment: Experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        db_models: List[WorkflowSchemaDbModel] = WorkflowSchemaDbModel.objects(experiment=experiment)

        if not db_models:
            raise NoLabsException(ErrorCodes.workflow_not_found)

        workflow_schema_model = db_models[0].get_workflow_value()

        workflow = Workflow.create_from_schema(
            workflow_schema_model=workflow_schema_model,
            functions=functions)

        await workflow.execute(terminate=True)
