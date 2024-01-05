import os
from typing import Dict

import numpy
import torch
import torch.nn as nn
import numpy as np
from argparse import Namespace
from esm import FastaBatchedDataset, pretrained  # type: ignore

import requests

from solubility.api_models import RunSolubilityPredictionRequest, RunSolubilityPredictionResponse
from solubility.loggers import Log
from solubility.settings import Settings

current_dir = os.path.dirname(os.path.abspath(__file__))
models_dir = os.path.join(current_dir, 'models')


Log.cuda_available(torch.cuda.is_available())


class SimpleSolubilityMultiLayerPerceptron(nn.Module):
    # Current MLP is build on top of ESM-2 150M
    def __init__(self, input_size = 640, hidden_sizes = [320, 160, 80, 40], output_size = 1, dropout_rate = 0.2):
        super(SimpleSolubilityMultiLayerPerceptron, self).__init__()
        layers = []
        prev_size = input_size
        for size in hidden_sizes:
            layers.append(nn.Linear(prev_size, size))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(dropout_rate))  # Add dropout layer
            prev_size = size
        self.hidden_layers = nn.Sequential(*layers)
        self.output_layer = nn.Linear(prev_size, output_size)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = self.hidden_layers(x)
        x = self.output_layer(x)
        x = self.sigmoid(x)
        return x


class ESM2EmbeddingGenerator:
    def __init__(self, model_name):
        gpu = torch.cuda.is_available()
        self.model_name = model_name
        self.model = None
        self.tokenizer = None
        self.gpu = gpu
        if gpu:
            self.device = torch.device('cuda')
        else:
            self.device = torch.device('cpu')

    def load_model(self):
        self.model, self.alphabet = pretrained.load_model_and_alphabet(self.model_name)
        self.model.to(self.device)
        self.model.eval()

    def predict(self, protein_sequence: str):
        include = "mean"
        repr_layers = [-1]
        truncation_seq_length = 1022
        nogpu = False
        toks_per_batch = 4096

        args = Namespace(
            repr_layers=repr_layers,
            include=include,
            truncation_seq_length=truncation_seq_length,
            nogpu=nogpu,
            toks_per_batch=toks_per_batch
        )

        if self.device == torch.device('cuda'):
            self.model = self.model.cuda()
            Log.transferred_models_to_gpu()

        dataset = FastaBatchedDataset(["mockId"], [protein_sequence])
        batches = dataset.get_batch_indices(args.toks_per_batch, extra_toks_per_seq=1)
        data_loader = torch.utils.data.DataLoader(
            dataset, collate_fn=self.alphabet.get_batch_converter(args.truncation_seq_length), batch_sampler=batches
        )

        return_contacts = "contacts" in args.include

        repr_layers = [(i + self.model.num_layers + 1) % (self.model.num_layers + 1) for i in args.repr_layers]

        with torch.no_grad():
            for batch_idx, (labels, strs, toks) in enumerate(data_loader):

                if self.device == torch.device('cuda'):
                    toks = toks.to(device="cuda", non_blocking=True)

                out = self.model(toks, repr_layers=repr_layers, return_contacts=return_contacts)

                logits = out["logits"].to(device="cpu")
                representations = {
                    layer: t.to(device="cpu") for layer, t in out["representations"].items()
                }

                for i, label in enumerate(labels):
                    label = label.split()[0]
                    result = {"label": label}
                    truncate_len = min(args.truncation_seq_length, len(strs[i]))
                    result["mean_representations"] = {
                        layer: t[i, 1: truncate_len + 1].mean(0).clone() \
                        for layer, t in representations.items()
                    }
                    return result["mean_representations"][30]


class SolubilityPrediction:
    def __init__(self, model_name, gpu):
        self.model_name = model_name
        self.model = None
        self.tokenizer = None
        self.gpu = gpu
        if gpu:
            self.device = torch.device('cuda')
        else:
            self.device = torch.device('cpu')
        self.embedding_model = None

    def load_model(self):
        Log.loading_solubility_model()
        self.model = SimpleSolubilityMultiLayerPerceptron().to(self.device)
        self._load_model_state_dict()
        self.model.eval()
        Log.solubility_has_been_loaded()

    def set_embedding_model(self, embedding_model):
        self.embedding_model = embedding_model

    def _load_model_state_dict(self):
        # URL of the model's .pth file on Hugging Face Model Hub
        model_url = "https://huggingface.co/thomasshelby/solubility_model/resolve/main/solubility_model.pth"

        # Path where you want to save the downloaded .pth file
        local_path = os.path.join(models_dir, "solubility_model.pth")

        if not os.path.exists(local_path):
            with open(local_path, 'wb') as f:
                response = requests.get(model_url)
                f.write(response.content)

        self.model.load_state_dict(torch.load(local_path, map_location=self.device))

    def predict(self, sequence: str):
        Log.making_solubility_predictions()
        embedding = self.embedding_model.predict(sequence)
        embedding = embedding.to(self.device)
        outputs = self.model(embedding.float())
        Log.successfully_predicted_solubility()
        return outputs.item()


def run_solubility_predictions(request: RunSolubilityPredictionRequest) -> RunSolubilityPredictionResponse:
    try:
        model = SolubilityPrediction(model_name='solubility', gpu=torch.cuda.is_available())
        model.load_model()
        if not model.embedding_model:
            settings = Settings()
            emb_model = ESM2EmbeddingGenerator(settings.model_weights_name)
            emb_model.load_model()
            model.set_embedding_model(emb_model)

        result = model.predict(request.proteinSequence)
        return RunSolubilityPredictionResponse(solubleConfidence=result, errors=[])
    except Exception:
        Log.exception()
        return RunSolubilityPredictionResponse(solubleConfidence=None,
                                                 errors=['Internal error in solubility prediction occured'])
