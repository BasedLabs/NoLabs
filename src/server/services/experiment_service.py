from abc import abstractmethod
from typing import Dict, Type, Union, Optional, List
import uuid

import os
import datetime
import json

from src.server.services.inference_service import create_pipeline, create_model, get_models_from_config
from src.server.services.savers import FileSaverFactory
from src.server.services.loaders import FileLoaderFactory, load_experiment_names


dirname = os.path.dirname
# Base directory for storing experiment results
EXPERIMENTS_DIR = dirname(dirname(dirname(os.path.abspath(__file__))) + "/experiments/")
PROTEIN_EXPERIMENTS_DIR = EXPERIMENTS_DIR + "/proteins"
DTI_EXPERIMENTS_DIR = EXPERIMENTS_DIR + "/drug_discovery"

# Helper function to ensure the base directory exists
def ensure_base_directory():
    if not os.path.exists(EXPERIMENTS_DIR):
        os.makedirs(EXPERIMENTS_DIR)
    if not os.path.exists(PROTEIN_EXPERIMENTS_DIR):
        os.makedirs(PROTEIN_EXPERIMENTS_DIR)
    if not os.path.exists(DTI_EXPERIMENTS_DIR):
        os.makedirs(DTI_EXPERIMENTS_DIR)

class BaseExperiment:

    @abstractmethod
    def run(self, *args, **kwargs):
        pass

class ProteinPropertyPrediction(BaseExperiment):

    def __init__(self, use_gpu=False, is_test=True):
        self.task2results_map = {
            "folding": "folding_prediction.pdb",
            "gene_ontology": "gene_ontology.json",
            "localisation": "cell_localisation.json",
            "solubility": "metadata"
        }
        self.pipeline = create_pipeline(use_gpu, is_test)
        models_metadata = get_models_from_config(is_test)
        for model_metadata in models_metadata:
            if model_metadata["task"] in self.task2results_map.keys():
                model = create_model(model_metadata, use_gpu)
                self.pipeline.add_model(model)

    def run(self, sequence: str, experiment_id: str):
        if not experiment_id:
            experiment_id = str(uuid.uuid4())
        for task, result_file in self.task2results_map:
            self.run_by_task(task, sequence, experiment_id, result_file)
        return experiment_id

    def run_by_task(self, task: str, sequence: str, experiment_id: str, result_file: str):
        model = self.pipeline.get_model_by_task(task)
        result = model.predict(sequence)
        self.store_result(experiment_id, result, result_file)

    def store_result(self, experiment_id: str, result, filename: str):
        file_saver = FileSaverFactory().get_saver(filename)
        experiment_dir = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id)
        file_saver.save(result, experiment_dir, filename)

    def save_experiment_metadata(experiment_id: str, experiment_name: str):
        metadata = {
            "id": experiment_id,
            "name": experiment_name,
            "date": datetime.now().isoformat()
        }
        
        metadata_path = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id, "metadata.json")
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=4)


    @classmethod
    def load_result(self, experiment_id: str, model_task: str):
        filename = self.task2results_map[model_task]
        experiment_dir = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id)
        loader = FileLoaderFactory().get_loader(filename)
        loaded_content = loader.load(experiment_dir, filename)
        return loaded_content

    @classmethod
    def load_experiment_names():
        return load_experiment_names(PROTEIN_EXPERIMENTS_DIR)


class DrugDiscovery(BaseExperiment):

    def __init__(self, use_gpu, is_test):
        # "skip" file fornat is used for skipping since the data saving procedure is written somewhere else
        # Here it's written in dti model predict method
        self.task2results_map = {
            "dti": "skip",
        }
        self.pipeline = create_pipeline(use_gpu, is_test)
        models_metadata = get_models_from_config(is_test)
        for model_metadata in models_metadata:
            if model_metadata["task"] in self.task2results_map.keys():
                model = create_model(model_metadata, use_gpu)
                self.pipeline.add_model(model)
    
    def run(self, ligand_files: List[str], protein_files: str, experiment_id):
        if not experiment_id:
            experiment_id = str(uuid.uuid4())
        model = self.pipeline.get_model_by_task("dti")
        experiment_dir = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id)
        model.predict(ligand_files, protein_files, experiment_dir)
        return experiment_id

    @classmethod
    def load_result(self, experiment_id: str, model_task: str):
        filename = self.task2results_map[model_task]
        experiment_dir = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id)
        loader = FileLoaderFactory().get_loader(filename)
        loaded_content = loader.load(experiment_dir, filename)
        return loaded_content

    @classmethod
    def load_result(self, experiment_id: str, model_task: str):
        filename = self.task2results_map[model_task]
        experiment_dir = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id)
        loader = FileLoaderFactory().get_loader(filename)
        loaded_content = loader.load(experiment_dir, filename)
        return loaded_content

    @classmethod
    def load_experiment_names():
        return load_experiment_names(DTI_EXPERIMENTS_DIR)
    
    def save_experiment_metadata(experiment_id: str, experiment_name: str):
        metadata = {
            "id": experiment_id,
            "name": experiment_name,
            "date": datetime.now().isoformat()
        }
        
        metadata_path = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id, "metadata.json")
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=4)

