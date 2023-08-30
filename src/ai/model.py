import torch
from torch import Tensor
from transformers import AutoModelForSequenceClassification, AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37

import numpy as np
import requests

from src.ai.exceptions.model_not_loaded_ex import ModelNotLoadedException
from src.ai.custom_models.custom_models import SimpleGOMultiLayerPerceptron, SimpleSolubilityMultiLayerPerceptron

from esm import Alphabet, FastaBatchedDataset, ProteinBertModel, pretrained, MSATransformer

from typing import List
import os

dirname = os.path.dirname

class BaseModel:
    def __init__(self, model_name: str, gpu: bool, model_task = ""):
        self.model_name = model_name
        self.model_task = model_task
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
    def __init__(self, model_name: str, gpu: bool, model_task = ""):
        super().__init__(model_name, gpu, model_task)

    def set_labels(self, labels: List[str]):
        self.labels = labels

    def load_model(self):
        self.model = AutoModelForSequenceClassification.from_pretrained(self.model_name)
        self.tokenizer = AutoTokenizer.from_pretrained(self.model_name)
        return super().load_model()

    def _raw_inference(self, sequence: str) -> Tensor:
        inputs = self.tokenizer(sequence, return_tensors='pt', padding=True, truncation=True)
        with torch.no_grad():
            outputs = self.model(**inputs)
            probabilities = torch.nn.functional.softmax(outputs.logits, dim=-1).tolist()[0]
        return probabilities

    def predict(self, sequence: str) -> List[List[tuple[str, int]]]:
        if not self.tokenizer or not self.model:
            raise ModelNotLoadedException()

        probabilities = self._raw_inference(sequence)
        prob_table = list(zip(self.labels, probabilities))
        return prob_table

class Folding(BaseModel):

    def __init__(self, model_name: str, gpu: bool, model_task = ""):
        super().__init__(model_name, gpu, model_task)
        self.model_name = model_name
        self.gpu = gpu

    def load_model(self):
        self.tokenizer = AutoTokenizer.from_pretrained(self.model_name)
        self.model = EsmForProteinFolding.from_pretrained(self.model_name, low_cpu_mem_usage=True)
        if self.gpu:
            self.model = self.model.cuda()

        return super().load_model()

    def convert_outputs_to_pdb(self, outputs) -> List[str]:
        final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)
        outputs = {k: v.to("cpu").numpy() for k, v in outputs.items()}
        final_atom_positions = final_atom_positions.cpu().numpy()
        final_atom_mask = outputs["atom37_atom_exists"]
        pdbs = []
        for i in range(outputs["aatype"].shape[0]):
            aa = outputs["aatype"][i]
            pred_pos = final_atom_positions[i]
            mask = final_atom_mask[i]
            resid = outputs["residue_index"][i] + 1
            pred = OFProtein(
                aatype=aa,
                atom_positions=pred_pos,
                atom_mask=mask,
                residue_index=resid,
                b_factors=outputs["plddt"][i],
                chain_index=outputs["chain_index"][i] if "chain_index" in outputs else None,
            )
            pdbs.append(to_pdb(pred))
        return pdbs


    def _raw_inference(self, sequence: str):

        tokenized_input = self.tokenizer([sequence], return_tensors="pt", add_special_tokens=False)['input_ids']
        
        if self.gpu:
            tokenized_input = tokenized_input.cuda()

        output = None

        with torch.no_grad():
            output = self.model(tokenized_input)

        return output

    def predict(self, sequence: str) -> List[str]:
        if not self.tokenizer or not self.model:
            raise ModelNotLoadedException()
        output = self._raw_inference(sequence)

        pdbs = self.convert_outputs_to_pdb(output)

        return "".join(pdbs)


class ESM2EmbeddingGenerator(BaseModel):
    def __init__(self, model_name, gpu, model_task = ""):
        super().__init__(model_name, gpu)
        if gpu:
            self.device = torch.device('cuda')
        else:
            self.device = torch.device('cpu')

    def load_model(self):
        self.model, self.alphabet = pretrained.load_model_and_alphabet(self.model_name)
        self.model.eval()

    def predict(self, protein_sequence: str):
        toks = self.alphabet.encode(protein_sequence)
        if self.device == torch.device('cuda'):
            toks = toks.to(device="cuda", non_blocking=True)
            self.model.to(self.device)
        repr_layers = [(i + self.model.num_layers + 1) % (self.model.num_layers + 1) for i in [-1]]
        out = self.model(toks, repr_layers=repr_layers, return_contacts=True)
        logits = out["logits"].to(device="cpu")
        representations = {
            layer: t.to(device="cpu") for layer, t in out["representations"].items()
        }
        result = {}
        result["mean_representations"] = {
                layer: t[i, 1 : truncate_len + 1].mean(0).clone()
                for layer, t in representations.items()
            }
        # For now it's 30 cus we using 150M ESM2
        # TODO: refactor on the weekend
        return result["mean_representations"][30]


class GeneOntologyPrediction(BaseModel):
    """
    Look for a better model there https://www.kaggle.com/competitions/cafa-5-protein-function-prediction
    """

    def __init__(self, model_name, gpu, model_task = ""):
        super().__init__(model_name, gpu)
        if gpu:
            self.device = torch.device('cuda')
        else:
            self.device = torch.device('cpu')
        self.embedding_model = None

    def load_model(self):
        self.model = SimpleGOMultiLayerPerceptron(input_dim=640, num_classes=200).to(self.device)
        self._load_model_state_dict()
        self.labels = np.load(dirname(os.path.abspath(__file__)) + \
             "/custom_models/models/gene_ontology/go_labels_200.npy", allow_pickle=True)

    def set_embedding_model(self, embedding_model):
        self.embedding_model = embedding_model

    def _load_model_state_dict(self):
        # URL of the model's .pth file on Hugging Face Model Hub
        model_url = "https://huggingface.co/thomasshelby/go_prediction/resolve/main/go_model_150M.pth"

        # Path where you want to save the downloaded .pth file
        local_path = dirname(os.path.abspath(__file__)) + "/custom_models/models/gene_ontology/go_model_150M.pth"

        response = requests.get(model_url)
        with open(local_path, 'wb') as f:
            f.write(response.content)

        self.model.load_state_dict(torch.load(local_path, map_location=self.device))

    def predict(self, protein_sequence: str):
        embedding = self.embedding_model.predict(protein_id)
        raw_outputs = torch.nn.functional.sigmoid(self.model(embedding)).squeeze().detach().cpu().numpy()
        
        outputs = {label: confidence for label, confidence in zip(self.labels, raw_outputs)}
        return outputs


class SolubilityPrediction(BaseModel):
    def __init__(self, model_name, gpu, model_task = ""):
        super().__init__(model_name, gpu)
        if gpu:
            self.device = torch.device('cuda')
        else:
            self.device = torch.device('cpu')
        self.embedding_model = None

    def load_model(self):
        self.model = SimpleSolubilityMultiLayerPerceptron().to(self.device)
        self._load_model_state_dict()

    def set_embedding_model(self, embedding_model):
        self.embedding_model = embedding_model

    def _load_model_state_dict(self):
        # URL of the model's .pth file on Hugging Face Model Hub
        model_url = "https://huggingface.co/thomasshelby/solubility_model/resolve/main/solubility_model.pth"

        # Path where you want to save the downloaded .pth file
        local_path = dirname(os.path.abspath(__file__)) + "/custom_models/models/solubility/solubility_model.pth"

        response = requests.get(model_url)
        with open(local_path, 'wb') as f:
            f.write(response.content)

        self.model.load_state_dict(torch.load(local_path, map_location=self.device))

    def predict(self, protein_id: str):
        embedding = self.embedding_model.predict(protein_id)
        embedding = embedding.to(self.device)
        outputs = model(embedding.float())
        predicted_labels = (outputs >= 0.5).long()

        return predicted_labels
