import datetime
import logging
from typing import List

from nolabs.exceptions import NoLabsException
from nolabs.workflow.component import Component
from nolabs.workflow.models import ComponentDbModel, JobErrorDbModel, InputPropertyErrorDbModel


class WorkflowExecutor:
    def __init__(self, logger: logging.Logger):
        self.logger = logger

    async def execute_single(self, component: Component):
        component_db_model: ComponentDbModel = ComponentDbModel.objects.with_id(component.id)

        input_changed = component.set_input_from_previous()

        jobs_setup_success = False
        last_jobs_count = component_db_model.last_jobs_count

        if input_changed or not component.jobs:
            component_db_model.input_parameter_dict = component.input_parameter_dict
            component_db_model.save()

            validation_errors = component.input_errors()

            if validation_errors:
                component_db_model.input_property_errors = [
                    InputPropertyErrorDbModel.create(
                        loc=e.loc,
                        msg=e.msg
                    ) for e in validation_errors
                ]
            else:
                try:
                    await component.setup_jobs()
                    component_db_model.jobs = component.jobs
                    component_db_model.last_jobs_count = len(component.jobs)
                    component_db_model.save()
                except Exception as e:
                    self.logger.exception('Exception occured in setup job')
                    if isinstance(e, NoLabsException):
                        component_db_model.last_exceptions = e.messages
                    else:
                        component_db_model.last_exceptions = [str(e)]
                    component_db_model.save()

                jobs_validation_errors = await component.prevalidate_jobs()

                if jobs_validation_errors:
                    component_db_model.jobs_errors = [JobErrorDbModel(
                        job_id=ve.job_id,
                        msg=ve.msg
                    ) for ve in jobs_validation_errors]
                else:
                    jobs_setup_success = True
                    component_db_model.jobs_errors = []

                component_db_model.save()

        if component.output_errors() or jobs_setup_success or any(
                [j for j in component.jobs if
                 j.updated_at >= component_db_model.last_executed_at]) or component_db_model.last_jobs_count != last_jobs_count:
            try:
                await component.execute()
                component_db_model.output_parameter_dict = component.output_parameter_dict
                component_db_model.last_executed_at = datetime.datetime.utcnow()
                component_db_model.save()
            except Exception as e:
                self.logger.exception('Exception occurred in execute job')
                if isinstance(e, NoLabsException):
                    component_db_model.last_exceptions = e.messages
                else:
                    component_db_model.last_exceptions = [str(e)]

                component_db_model.save()

        component_db_model.last_execute_try_at = datetime.datetime.utcnow()
        component_db_model.save()

    async def execute(self, components: List[Component]):
        async def execute(component: Component):
            component_db_model: ComponentDbModel = ComponentDbModel.objects.with_id(component.id)

            for previous in component.previous:
                await execute(previous)

                validation_errors = previous.output_errors()

                if validation_errors:
                    component_db_model.input_property_errors = [InputPropertyErrorDbModel(loc=err.loc, msg=err.msg) for
                                                                err in validation_errors]
                    component_db_model.save()
                    return

            input_changed = component.set_input_from_previous()

            jobs_setup_success = False
            last_jobs_count = component_db_model.last_jobs_count

            if input_changed or not component.jobs:
                component_db_model.input_parameter_dict = component.input_parameter_dict
                component_db_model.save()

                validation_errors = component.input_errors()

                if validation_errors:
                    component_db_model.input_property_errors = [
                        InputPropertyErrorDbModel.create(
                            loc=e.loc,
                            msg=e.msg
                        ) for e in validation_errors
                    ]
                else:
                    try:
                        await component.setup_jobs()
                        component_db_model.jobs = component.jobs
                        component_db_model.last_jobs_count = len(component.jobs)
                        component_db_model.save()
                    except Exception as e:
                        self.logger.exception('Exception occured in setup job')
                        if isinstance(e, NoLabsException):
                            component_db_model.last_exceptions = e.messages
                        else:
                            component_db_model.last_exceptions = [str(e)]
                        component_db_model.save()

                    jobs_validation_errors = await component.prevalidate_jobs()

                    if jobs_validation_errors:
                        component_db_model.jobs_errors = [JobErrorDbModel(
                            job_id=ve.job_id,
                            msg=ve.msg
                        ) for ve in jobs_validation_errors]
                    else:
                        jobs_setup_success = True
                        component_db_model.jobs_errors = []

                    component_db_model.save()

            if component.output_errors() or jobs_setup_success or any(
                    [j for j in component.jobs if
                     j.updated_at >= component_db_model.last_execute_try_at]) or component_db_model.last_jobs_count != last_jobs_count:
                try:
                    await component.execute()
                    component_db_model.output_parameter_dict = component.output_parameter_dict
                    component_db_model.last_executed_at = datetime.datetime.utcnow()
                    component_db_model.save()
                except Exception as e:
                    self.logger.exception('Exception occurred in execute job')
                    if isinstance(e, NoLabsException):
                        component_db_model.last_exceptions = e.messages
                    else:
                        component_db_model.last_exceptions = [str(e)]

                    component_db_model.save()

            component_db_model.last_execute_try_at = datetime.datetime.utcnow()
            component_db_model.save()

        for component in components:
            await execute(component=component)
