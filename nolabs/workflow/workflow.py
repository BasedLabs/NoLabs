from dataclasses import dataclass
from typing import List

from nolabs.refined.domain.models.common import Experiment
from nolabs.workflow.application.models import WorkflowSchemaDbModel
from nolabs.workflow.component import PythonComponent
from nolabs.workflow.workflow_schema import WorkflowSchemaModel, JobValidationError


@dataclass
class WorkflowValidationError:
    msg: str


class Workflow:
    async def execute(self, workflow_schema: WorkflowSchemaModel, components: List[PythonComponent]):
        experiment = Experiment.objects.with_id(workflow_schema.experiment_id)
        db_model: WorkflowSchemaDbModel = WorkflowSchemaDbModel.objects(experiment=experiment).first()

        async def execute(component: PythonComponent):
            for previous in component.previous:
                await execute(previous)

                validation_errors = previous.validate_output()

                if validation_errors:
                    workflow_schema.get_wf_component(previous.id).error = ', '.join(
                        [ve.msg for ve in validation_errors])
                    workflow_schema.valid = False
                    db_model.set_workflow_value(workflow_schema)
                    db_model.save()

                return

            input_changed = component.set_input_from_previous()

            if not input_changed:
                return

            validation_errors = component.validate_input()

            if validation_errors:
                workflow_schema.get_wf_component(component.id).error = ', '.join(
                    [ve.msg for ve in validation_errors])
                db_model.set_workflow_value(workflow_schema)
                db_model.save()
                return

            await component.setup_jobs()

            jobs_validation_errors = await component.prevalidate_jobs()

            if jobs_validation_errors:
                workflow_schema.get_wf_component(component.id).error = [JobValidationError(
                    job_id=ve.job_id,
                    msg=ve.msg
                ) for ve in jobs_validation_errors]
            else:
                workflow_schema.get_wf_component(component.id).jobs_errors = []
                db_model.set_workflow_value(workflow_schema)
                db_model.save()

            await component.execute()

        for component in components:
            await execute(component=component)
