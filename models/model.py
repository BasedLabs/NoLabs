import torch
from transformers import AutoModelForSequenceClassification, AutoTokenizer

class BaseModel:
    def __init__(self, model_name):
        self.model_name = model_name
        self.model = None
        self.tokenizer = None

    def load_model(self):
        self.model = AutoModelForSequenceClassification.from_pretrained(self.model_name)
        self.tokenizer = AutoTokenizer.from_pretrained(self.model_name)

    #Method to get raw model outputs
    def raw_inference(self, input):
        pass

    # Method to return raw outputs in the desired format
    def predict(self, input):
        pass

class LocalisationModel(BaseModel):
    def __init__(self, model_name):
        super().__init__(model_name)
    
    def set_labels(self, labels):
        self.labels = labels

    def raw_inference(self, sequence):
        inputs = self.tokenizer(sequence, return_tensors='pt', padding=True, truncation=True)
        with torch.no_grad():
            outputs = self.model(**inputs)
            probabilities = torch.nn.functional.softmax(outputs.logits, dim=-1).tolist()[0]
        return probabilities

    def predict(self, sequence):
        probabilities = self.raw_inference(sequence)
        prob_table = list(zip(self.labels, probabilities))

        return prob_table

        

class FunctionPrediction(BaseModel):
    def __init__(self, model_name):
        super().__init__(model_name)

    def predict(self):
        #TODO: complete function prediction from https://github.com/kexinhuang12345/DeepPurpose
        pass