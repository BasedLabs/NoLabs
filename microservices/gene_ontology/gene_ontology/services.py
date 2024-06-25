import os
from typing import Dict

import numpy
import torch
import torch.nn as nn
import numpy as np
from argparse import Namespace
from esm import FastaBatchedDataset, pretrained  # type: ignore

import requests

from gene_ontology.api_models import RunGeneOntologyPredictionRequest, RunGeneOntologyPredictionResponse, \
    GoConfidenceResponse
from gene_ontology.loggers import logger
from gene_ontology.settings import Settings
from gene_ontology.shared import memory_manager, JobStatusEnum

current_dir = os.path.dirname(os.path.abspath(__file__))
models_dir = os.path.join(current_dir, 'models')


logger.cuda_available(torch.cuda.is_available())


class SimpleGOMultiLayerPerceptron(nn.Module):
    def __init__(self, input_dim, num_classes):
        super(SimpleGOMultiLayerPerceptron, self).__init__()
        self.linear1 = nn.Linear(input_dim, 500)
        self.activation1 = nn.ReLU()
        self.linear2 = nn.Linear(500, 300)
        self.activation2 = nn.ReLU()
        self.linear3 = nn.Linear(300, num_classes)

    def forward(self, x):
        x = self.linear1(x)
        x = self.activation1(x)
        x = self.linear2(x)
        x = self.activation2(x)
        x = self.linear3(x)
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
            logger.transferred_models_to_gpu()

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


class GeneOntologyPrediction:
    """
    Look for a better model there https://www.kaggle.com/competitions/cafa-5-protein-function-prediction
    """

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
        logger.loading_gene_ontology_model()
        self.model = SimpleGOMultiLayerPerceptron(input_dim=640, num_classes=200).to(self.device)
        self._load_model_state_dict()
        self.model.eval()
        self.labels = np.load(os.path.join(models_dir, "go_labels_200.npy"), allow_pickle=True)
        logger.gene_ontology_has_been_loaded()

    def set_embedding_model(self, embedding_model):
        self.embedding_model = embedding_model

    def _load_model_state_dict(self):
        # URL of the model's .pth file on Hugging Face Model Hub
        model_url = "https://huggingface.co/thomasshelby/go_prediction/resolve/main/go_model_150M.pth"

        # Path where you want to save the downloaded .pth file
        local_path = os.path.join(models_dir, "go_model_150M.pth")

        if not os.path.exists(local_path):
            with open(local_path, 'wb') as f:
                response = requests.get(model_url)
                f.write(response.content)

        self.model.load_state_dict(torch.load(local_path, map_location=self.device))

    def predict(self, protein_sequence: str) -> Dict[str, numpy.float32]:
        logger.making_gene_ontology_predictions()
        embedding = self.embedding_model.predict(protein_sequence)
        embedding = embedding.to(self.device)
        # type: ignore
        raw_outputs = torch.nn.functional.sigmoid(self.model(embedding)).squeeze().detach().cpu().numpy()

        outputs: Dict[str, numpy.float32] = {label: confidence for label, confidence in zip(self.labels, raw_outputs)}
        logger.successfully_predicted_gene_ontology()
        return outputs


def run_gene_ontology_prediction(request: RunGeneOntologyPredictionRequest) -> RunGeneOntologyPredictionResponse:
    try:
        memory_manager.change_status(str(request.job_id), JobStatusEnum.running)
        model = GeneOntologyPrediction(model_name='gene_ontology',
                                       gpu=torch.cuda.is_available())
        model.load_model()
        if not model.embedding_model:
            settings = Settings()
            emb_model = ESM2EmbeddingGenerator(settings.model_weights_name)
            emb_model.load_model()
            model.set_embedding_model(emb_model)

        result = model.predict(request.amino_acid_sequence)
        native_float_result = {}
        for key in result.keys():
            native_float_result[key] = result[key].item()
        return RunGeneOntologyPredictionResponse(go_confidence=[
            GoConfidenceResponse(name=name, confidence=confidence) for name, confidence in result.items()
        ], errors=[])
    except Exception:
        logger.exception()
        return RunGeneOntologyPredictionResponse(go_confidence=[],
                                                 errors=['Internal error in gene ontology prediction occured'])
    finally:
        memory_manager.change_status(str(request.job_id), JobStatusEnum.idle)
