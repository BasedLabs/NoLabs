import uuid
from typing import List

from airflow import DAG

from nolabs.application.workflow import ComponentBag


def create_dynamic_dag(workflow_id: uuid.UUID, components: List[ComponentBag]):
    dag = DAG(
        dag_id=str(workflow_id),
        default_args={"owner": "airflow"},
        schedule_interval=None
    )

    tasks = {}

    for component in components:
        # Assuming each component has a unique name
        if component.JobsSetupOperatorType:
            setup_task_name = f"{component.component_name}_{component.JobsSetupOperatorType.__name__}"

            setup_task = component.JobsSetupOperatorType(component_id=task_id=setup_task_name, dag=dag)
            tasks[setup_task_name] = setup_task

        job_task = component.JobOperatorType(task_id=f"{component.component_name}_job", dag=dag)
        tasks[f"{component.component_name}_job"] = job_task

        if component.OutputOperatorType:
            output_task = component.OutputOperatorType(task_id=f"{component.component_name}_output", dag=dag)
            tasks[f"{component.component_name}_output"] = output_task

    # Set dependencies
    for component in components:
        if component.JobsSetupOperatorType:
            setup_task = tasks[f"{component.component_name}_setup"]
            job_task = tasks[f"{component.component_name}_job"]
            setup_task >> job_task

        if component.OutputOperatorType:
            job_task = tasks[f"{component.component_name}_job"]
            output_task = tasks[f"{component.component_name}_output"]
            job_task >> output_task

        # Example for adding dependencies between different components, modify as needed
        # if component.depends_on:
        #     for dependency in component.depends_on:
        #         dependent_task = tasks[f"{dependency}_output"]
        #         setup_task = tasks[f"{component.component_name}_setup"]
        #         dependent_task >> setup_task

    return dag