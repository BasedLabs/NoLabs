from typing import List

import torch
from torch import Tensor

from localisation.api_models import RunLocalisationPredictionResponse, RunLocalisationPredictionRequest
from transformers import AutoModelForSequenceClassification, AutoTokenizer, EsmForProteinFolding
from localisation.settings import Settings


class ClassificationModel:
    def __init__(self):
        self.model_name = Settings().model_name
        self.model = AutoModelForSequenceClassification.from_pretrained(self.model_name)
        self.tokenizer = AutoTokenizer.from_pretrained(self.model_name)

    def load_model(self):
        self.model.eval()

    def _raw_inference(self, sequence: str) -> Tensor:
        inputs = self.tokenizer(sequence, return_tensors='pt', padding=True, truncation=True)
        with torch.no_grad():
            outputs = self.model(**inputs)
            probabilities = torch.nn.functional.softmax(outputs.logits, dim=-1).tolist()[0]
        return probabilities

    def predict(self, request: RunLocalisationPredictionRequest) -> RunLocalisationPredictionResponse:
        probabilities = self._raw_inference(request.amino_acid_sequence)
        response = RunLocalisationPredictionResponse(
            errors=[],
            cytosolic_proteins=probabilities[0],
            mitochondrial_proteins=probabilities[1],
            nuclear_proteins=probabilities[2],
            other_proteins=probabilities[3],
            extracellular_secreted_proteins=probabilities[4],
        )
        return response


def run_localisation(request: RunLocalisationPredictionRequest) -> RunLocalisationPredictionResponse:
    classification_model = ClassificationModel()
    return classification_model.predict(request)
