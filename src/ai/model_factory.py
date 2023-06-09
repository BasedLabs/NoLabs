import torch
from typing import Dict, Union

from .exceptions.unknown_model_ex import UnknownModelException
from .model import ClassificationModel, Folding, BaseModel


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

    raise UnknownModelException()
