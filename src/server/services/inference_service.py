import json
import os
from typing import Dict

from src.ai.model_factory import create_model
from src.ai.pipeline import Pipeline


def create_pipeline(use_gpu=False, is_test = False) -> Pipeline:
    pipeline = Pipeline()
    models_metadata = get_models_from_config(is_test)

    for model_metadata in models_metadata:
        model = create_model(model_metadata, use_gpu)
        pipeline.add_model(model)

    return pipeline


def get_localisation_output(pipeline, amino_acid_sequence: str) -> Dict:
    assert amino_acid_sequence
    model = pipeline.get_model_by_task("localisation")
    return model.predict(amino_acid_sequence)


def get_folding_output(pipeline, amino_acid_sequence: str) -> str:
    assert amino_acid_sequence
    model = pipeline.get_model_by_task("folding")
    return model.predict(amino_acid_sequence)

def get_gene_ontology_output(pipeline, protein_id: str) -> Dict:
    model = pipeline.get_model_by_task("gene_ontology")
    return model.predict(protein_id)

def get_pipeline_output(pipeline, amino_acid_sequence: str) -> str:
    assert amino_acid_sequence
    return pipeline.predict(amino_acid_sequence)


def read_json_file(file_path):
    with open(file_path, 'r') as json_file:
        data = json.load(json_file)
    return data


def read_config():
    dirname = os.path.dirname
    root_directory = dirname(dirname(dirname(dirname(os.path.abspath(__file__)))))
    file_path = os.path.join(root_directory, 'config.json')
    data = read_json_file(file_path)
    return data

def read_test_config():
    dirname = os.path.dirname
    root_directory = dirname(dirname(dirname(dirname(os.path.abspath(__file__)))))
    file_path = os.path.join(root_directory, 'test/ai/mock_resources/config.json')
    data = read_json_file(file_path)
    return data

def get_models_from_config(is_test):
    data = None
    if is_test:
        data = read_test_config()
    else:
        data = read_config()
    models = data.get('models', [])
    return models
