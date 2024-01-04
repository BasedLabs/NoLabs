import os
from typing import List
from nolabs.ai.model import BaseModel

import requests

class APIFolding(BaseModel):

    def __init__(self, model_name: str, gpu: bool, model_task = ""):
        super().__init__(model_name, gpu, model_task)
        self.model_name = model_name
        self.gpu = gpu

    def load_model(self):
        pass

    def read_pdb_file(self, file_path):
        try:
            with open(file_path, 'r') as file:
                pdb_contents = file.read()
                return "".join(pdb_contents.split())
        except FileNotFoundError:
            print(f"Error: File '{file_path}' not found.")
            return None

    def predict(self, sequence: str) -> List[str]:
        url = "https://api.esmatlas.com/foldSequence/v1/pdb/" 
        response = requests.post(url, data=sequence, verify=False)
        
        return response.text


