import torch
from typing import Dict, Union

from .exceptions.unknown_model_ex import UnknownModelException
from .model import ClassificationModel, Folding, SolubilityPrediction, \
    GeneOntologyPrediction, ESM2EmbeddingGenerator, DrugTargetInteraction, BaseModel
from test.ai.mock_model import FakeFolding


def create_model(model_metadata: Dict[str, str], use_gpu: bool = False) -> Union[BaseModel, None]:
    assert 'name' in model_metadata
    assert 'task' in model_metadata

    model_name = model_metadata.get('name')
    model_type = model_metadata.get('type')

    device = "cpu"

    if use_gpu:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    if device == "cuda":
        use_gpu = True

    if model_type == "classification":
        model_labels = model_metadata.get('labels', [])

        model_task = ""

        if 'task' in model_metadata:
            model_task = model_metadata["task"]

        model = ClassificationModel(model_name=model_name, gpu=use_gpu, model_task=model_task)
        model.load_model()
        model.set_labels(model_labels)

        return model

    if model_type == "folding":
        
        model_task = ""

        if 'task' in model_metadata:
            model_task = model_metadata["task"]

        model = Folding(model_name=model_name, gpu=use_gpu, model_task=model_task)
        model.load_model()
        
        return model

    # adding for testing purposes (folding is not fast :()), will substitute with an api call option
    if model_type == "fakefolding":

        model_task = ""

        if 'task' in model_metadata:
            model_task = model_metadata["task"]

        model = FakeFolding(model_name=model_name, gpu=use_gpu, model_task=model_task)
        model.load_model()

        return model

    if model_type == "gene_ontology":

        model_task = ""

        if 'task' in model_metadata:
            model_task = model_metadata["task"]

        model = GeneOntologyPrediction(model_name=model_name, gpu=use_gpu, model_task=model_task)
        model.load_model()
        if not model.embedding_model:
            emb_model = ESM2EmbeddingGenerator(model_metadata['embedding_model'], gpu=use_gpu)
            emb_model.load_model()
            model.set_embedding_model(emb_model)

        return model

    if model_type == "solubility":

        model_task = ""

        if 'task' in model_metadata:
            model_task = model_metadata["task"]

        model = SolubilityPrediction(model_name=model_name, gpu=use_gpu, model_task=model_task)
        model.load_model()
        if not model.embedding_model:
            emb_model = ESM2EmbeddingGenerator(model_metadata['embedding_model'], gpu=use_gpu)
            emb_model.load_model()
            model.set_embedding_model(emb_model)

        return model

    if model_type == "dti":

        model_task = ""

        if 'task' in model_metadata:
            model_task = model_metadata["task"]

        model = DrugTargetInteraction(model_name=model_name, gpu=use_gpu, model_task=model_task)
        model.load_model()

        return model

    raise UnknownModelException()
