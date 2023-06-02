import torch
from transformers import AutoModelForSequenceClassification, AutoTokenizer

from src.ai.exceptions.model_not_loaded_ex import ModelNotLoadedException


class BaseModel:
    def __init__(self, model_name, gpu):
        self.model_name = model_name
        self.model = None
        self.tokenizer = None
        self.gpu = gpu

    def load_model(self):
        """Load model and tokenizer here"""
        pass

    # Method to get raw model outputs
    def _raw_inference(self, input: str):
        pass

    # Method to return raw outputs in the desired format
    def predict(self, input: str):
        pass


class ClassificationModel(BaseModel):
    def __init__(self, model_name, gpu):
        super().__init__(model_name, gpu)

    def set_labels(self, labels):
        self.labels = labels

    def load_model(self):
        self.model = AutoModelForSequenceClassification.from_pretrained(self.model_name)
        self.tokenizer = AutoTokenizer.from_pretrained(self.model_name)
        return super().load_model()

    def _raw_inference(self, sequence: str):
        inputs = self.tokenizer(sequence, return_tensors='pt', padding=True, truncation=True)
        with torch.no_grad():
            outputs = self.model(**inputs)
            probabilities = torch.nn.functional.softmax(outputs.logits, dim=-1).tolist()[0]
        return probabilities

    def predict(self, sequence: str):
        if not self.tokenizer or not self.model:
            raise ModelNotLoadedException()

        probabilities = self._raw_inference(sequence)
        prob_table = list(zip(self.labels, probabilities))
        return prob_table


class FunctionPrediction(BaseModel):
    def __init__(self, model_name, gpu):
        super().__init__(model_name, gpu)

    def predict(self, sequence):
        # TODO: complete function prediction from https://github.com/kexinhuang12345/DeepPurpose
        pass
