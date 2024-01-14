import requests
from msa_light.api_models import RunMsaPredictionRequest, RunMsaPredictionResponse
from msa_light.loggers import Log

__all__ = ['predict_msa_service']


def predict_msa_service(parameters: RunMsaPredictionRequest) -> RunMsaPredictionResponse:
    try:
        files = {'sequence_file': ('filename', parameters.fasta_contents)}

        # Send the POST request with the file content
        response = requests.post(parameters.api_url, files=files)

        # Check if the request was successful
        if response.status_code == 200:
            return RunMsaPredictionResponse(msa_contents= response.json()['alignment'])
        return RunMsaPredictionResponse(msa_contents=None,
                                            errors=['Unable to run esmfold. The api server returned a status code: ', str(response.status_code)])
    except Exception as e:
        Log.exception()
        return RunMsaPredictionResponse(msa_contents=None,
                                            errors=['Unable to run esmfold. Internal error', str(e)])
