from abc import abstractmethod
from typing import Dict, Type, Union, Optional, List
import uuid

import os
import datetime
import json

from src.server.services.inference_service import create_pipeline, create_model, get_models_from_config
from src.server.services.savers import FileSaverFactory
from src.server.services.loaders import FileLoaderFactory


dirname = os.path.dirname
# Base directory for storing experiment results
EXPERIMENTS_DIR = dirname(dirname(dirname(os.path.abspath(__file__))) + "/experiments/")

# Helper function to ensure the base directory exists
def ensure_base_directory():
    if not os.path.exists(EXPERIMENTS_DIR):
        os.makedirs(EXPERIMENTS_DIR)


# Helper function to save experiment metadata
def save_experiment_metadata(experiment_id: str, experiment_name: str):
    metadata = {
        "id": experiment_id,
        "name": experiment_name,
        "date": datetime.now().isoformat()
    }
    
    metadata_path = os.path.join(EXPERIMENTS_DIR, experiment_id, "metadata.json")
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=4)


class BaseExperiment:

    @abstractmethod
    def run(self, *args, **kwargs):
        pass

class ProteinPropertyPrediction(BaseExperiment):

    def __init__(self, use_gpu, is_test):
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
        for task, result_file in self.task2results_map:
            self.run_by_task(task, sequence, experiment_id, result_file)

    def run_by_task(self, task: str, sequence: str, experiment_id: str, result_file: str):
        model = self.pipeline.get_model_by_task(task)
        result = model.predict(sequence)
        self.store_result(experiment_id, result, result_file)

    def store_result(self, experiment_id: str, result, filename: str):
        file_saver = FileSaverFactory().get_saver(filename)
        experiment_dir = os.path.join(EXPERIMENTS_DIR, experiment_id)
        file_saver.save(result, experiment_dir, filename)

    @classmethod
    def load_result(self, experiment_id: str, model_task: str):
        filename = self.task2results_map[model_task]
        experiment_dir = os.path.join(EXPERIMENTS_DIR, experiment_id)
        loader = FileLoaderFactory().get_loader(filename)
        loaded_content = loader.load(experiment_dir, filename)
        return loaded_content



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
    
    def run_dti(pipeline, ligand_files: List[str], protein_file: str):
        model = pipeline.get_model_by_task("dti")
        model.predict(ligand_files, protein_file)


class ExperimentFactory:

    _experiments: Dict[str, Type[BaseExperiment]] = {}

    @classmethod
    def register_experiment(cls, key: str, experiment: Type[BaseExperiment]):
        if key in cls._experiments:
            raise ValueError(f"Experiment {key} already registered!")
        cls._experiments[key] = experiment

    @classmethod
    def get_experiment(cls, key: str) -> Type[BaseExperiment]:
        experiment = cls._experiments.get(key)
        if not experiment:
            raise ValueError(f"Experiment {key} not found!")
        return experiment


# Experiment results manager
class ExperimentResultsManager:

    _results_store: Dict[str, Union[str, Dict]] = {}

    @classmethod
    def get_result(cls, experiment_id: str) -> Optional[Union[str, Dict]]:
        return cls._results_store.get(experiment_id)

class ExperimentAPI:

    @staticmethod
    def run_experiment(name_or_id: str, *args, **kwargs) -> str:
        ExperimentClass = ExperimentFactory.get_experiment(name_or_id)
        experiment_instance = ExperimentClass()
        experiment_id = str(uuid.uuid4())
        result = experiment_instance.run(*args, **kwargs)

        return experiment_id  # Return the unique experiment ID

    @staticmethod
    def get_experiment_result(experiment_id: str) -> Optional[Union[str, Dict]]:
        return ExperimentResultsManager.get_result(experiment_id)