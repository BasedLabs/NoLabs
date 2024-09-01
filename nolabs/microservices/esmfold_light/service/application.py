import requests

from api_models import InferenceInput, InferenceOutput


def inference(param: InferenceInput) -> InferenceOutput:
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    response = requests.post(url, data=param.fasta_sequence, verify=False)

    return InferenceOutput(pdb_content=response.text)