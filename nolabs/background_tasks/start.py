from dotenv import load_dotenv

load_dotenv("infrastructure/.env")

from prefect import serve

from nolabs.background_tasks.jobs import (
    cleanup_orhpan_flow_run_ids,
    cleanup_orhpan_task_run_ids,
)

if __name__ == "__main__":
    serve(
        cleanup_orhpan_task_run_ids.to_deployment(
            name="cleanup_orphan_task_run_ids", interval=10
        ),
        cleanup_orhpan_flow_run_ids.to_deployment(
            name="cleanup_orphan_flow_run_ids", interval=10
        ),
    )
