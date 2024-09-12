import requests

from api_models import InferenceInput, InferenceOutput
from log import logger


def inference(param: InferenceInput) -> InferenceOutput:
    logger.info('Starting inference', extra={
        'fasta_truncated': param.fasta_sequence[0:5]
    })
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    response = requests.post(url, data=param.fasta_sequence, verify=False)
    logger.info('Finished inference', extra={
        'pdb_truncated': response.text[0:5]
    })

    return InferenceOutput(pdb_content=response.text)