from typing import Any

from airflow.utils.context import Context

from nolabs.application.workflow import ExecuteJobOperator




class RunEsmFoldLightOperator(ExecuteJobOperator):
    async def execute_async(self, context: Context) -> Any:
        import requests
        from nolabs.job_services.esmfold_light.operators.api_models import PredictFoldingJobRequest, PredictFoldingJobResponse

        input = PredictFoldingJobRequest(**self.get_input())

        url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
        response = requests.post(url, data=input.protein_sequence, verify=False)

        self.set_output(PredictFoldingJobResponse(pdb_content=response.text))
