import os
import csv
from datetime import date
from typing import List

import time

from .model import BaseModel

from src.server.services.experiments_structure_loader import ProteinLabExperimentsLoader
from src.server.services.progress import ProgressTracker

class Pipeline:
    def __init__(self, models: List[BaseModel] = [], progress_file='progress.json'):
        """
        models: pass the list of loaded models which inherit BaseModel class
        """
        print("pipeline receives models: ", models)
        self.models = models

    def get_model_names(self) -> List[str]:
        return [model.model_name for model in self.models]
    
    def get_model_tasks(self) -> List[str]:
        return [model.model_task for model in self.models]

    def get_model_by_task(self, model_task: str):
        """
        Gets the first model of the required type from the pipeline (TODO: get all the models of correct type)
        """
        for model in self.models:
            if model.model_task == model_task:
                return model

    def add_model(self, model: BaseModel):
        """
        model: loaded model which inherits the BaseModel class
        """
        self.models.append(model)

    def predict_single(self, sequence: str):
        """
        Make predictions for one sequence
        """
        output = {}
        for model in self.models:
            print(model.model_task)
            output[model.model_task] = model.predict(sequence)
        return output

    def predict_fasta(self, sequences: List[str],
            protein_ids: List[str],
            loader: ProteinLabExperimentsLoader,
            experiment_dir: str,):
        """
        Make predictions for multiple sequences
        """

        experiment_progress_tracker = ProgressTracker(experiment_dir, protein_ids)
        for sequence, protein_id in zip(sequences, protein_ids):
            output = self.predict_single(sequence)
            protein_dir = os.path.join(experiment_dir, protein_id)

            protein_progress_tracker = ProgressTracker(protein_dir, tasks=self.get_model_tasks())
            for task, result in output.items():
                loader.save_outputs(result=result, task=task, save_dir=protein_dir)
                protein_progress_tracker.update_progress(task)
            experiment_progress_tracker.update_progress(protein_id)

    def run(self, fasta_files_paths, loader: ProteinLabExperimentsLoader, experiment_dir=""):
        protein_ids, sequences = loader.get_sequences(fasta_files_paths)

        self.predict_fasta(protein_ids=protein_ids, 
                           sequences=sequences, 
                           loader=loader, 
                           experiment_dir=experiment_dir)


def check_folder_exists(filename):
    directory = os.path.dirname(filename)

    if os.path.exists(directory):
        print(f"The folder for '{filename}' does not exist. Creating one")
        os.makedirs(directory)
