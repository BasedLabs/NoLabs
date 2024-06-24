import datetime
import logging
from typing import List, Optional

from nolabs.exceptions import NoLabsException
from nolabs.application.workflow.component import Component
from nolabs.application.workflow.models import ComponentDbModel, JobErrorDbModel, InputPropertyErrorDbModel


class WorkflowExecutor:
    def __init__(self, logger: logging.Logger):
        self.logger = logger

    async def execute(self, components: List[Component], specific: Optional[Component] = None):
        async def execute(component: Component):
            component_db_model: ComponentDbModel = ComponentDbModel.objects.with_id(component.id)
            try:
                if not specific:
                    for previous in component.previous:
                        await execute(previous)

                        validation_errors = previous.output_errors()

                        if validation_errors:
                            component_db_model.input_property_errors = [InputPropertyErrorDbModel(loc=err.loc, msg=err.msg)
                                                                        for
                                                                        err in validation_errors]
                            component_db_model.save()
                            return

                input_changed = component.set_input_from_previous()

                last_jobs_count = component_db_model.last_jobs_count

                setup_job_issue = False

                if input_changed or not component.jobs:
                    component_db_model.input_parameter_dict = component.input_parameter_dict
                    component_db_model.save()

                    validation_errors = component.input_errors()

                    if not validation_errors:
                        await component.setup_jobs()
                        component_db_model.jobs = component.jobs
                        component_db_model.last_jobs_count = len(component.jobs)
                        component_db_model.input_property_errors = []
                        component_db_model.save()
                    else:
                        component_db_model.input_property_errors = [
                            InputPropertyErrorDbModel.create(
                                loc=e.loc,
                                msg=e.msg
                            ) for e in validation_errors
                        ]

                jobs_validation_errors = await component.jobs_setup_errors()

                def jobs_updated():
                    return any([j for j in component.jobs if
                                not j.inputs_updated_at or
                                not component_db_model.last_executed_at or
                                j.inputs_updated_at >= component_db_model.last_executed_at])

                if ((jobs_updated() or component_db_model.last_jobs_count != last_jobs_count or input_changed)
                        and not jobs_validation_errors
                        and not setup_job_issue):
                    await component.execute()
                    component_db_model.output_parameter_dict = component.output_parameter_dict
                    component_db_model.last_executed_at = datetime.datetime.utcnow()
                    component_db_model.input_property_errors = []
                    component_db_model.last_exceptions = []
                    component_db_model.save()
            except Exception as e:
                self.logger.exception('Exception occurred in execute job')
                if isinstance(e, NoLabsException):
                    component_db_model.last_exceptions = e.messages
                else:
                    component_db_model.last_exceptions = [str(e)]
                component_db_model.save()

        if specific:
            await execute(component=specific)
            return

        for component in components:
            await execute(component=component)

