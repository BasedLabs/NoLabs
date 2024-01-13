import requests

class APIFolding:

    def __init__(self):
        pass

    def predict(self, sequence: str) -> str:
        url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
        response = requests.post(url, data=sequence, verify=False)

        return response.text

