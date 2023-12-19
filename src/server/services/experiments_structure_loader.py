import json
import os
import shutil
from abc import ABC, abstractmethod
from datetime import datetime
from typing import List, Dict
from pathlib import Path
from werkzeug.datastructures import FileStorage

from src.server.services.loaders import FileLoaderFactory, DTILoader, FastaFileLoader
from src.server.services.savers import FileSaverFactory, FastaFileSaver
from src.server.settings import PROTEIN_EXPERIMENTS_DIR, DTI_EXPERIMENTS_DIR, CONFORMATIONS_EXPERIMENTS_DIR
from src.server.services.progress import get_progress, get_protein_progress

def _load_experiments_ids_names(experiments_dir) -> Dict:
    # List all directories inside the dti experiments folder
    experiment_ids = [d for d in os.listdir(experiments_dir) if os.path.isdir(os.path.join(experiments_dir, d))]

    result_dict = {}

    for experiment_id in experiment_ids:
        experiment_dir = os.path.join(experiments_dir, experiment_id)
        metadata_path = os.path.join(experiment_dir, 'metadata.json')
    
        # Check if metadata.json exists in the directory
        if os.path.exists(metadata_path):
            with open(metadata_path, 'r') as f:
                data = json.load(f)
                # Assuming the JSON structure has a key 'experiment_name' with the name of the experiment.
                # Adjust if the structure is different.
                experiment_name = data.get('name', None)
                experiment_date = data.get('date', None)

                #Get progress if exists
                experiment_progress = 0
                try:
                    experiment_progress = get_progress(experiment_dir)['progress']
                except:
                    print('Progress does not exist for experiment: ', experiment_id)

                if experiment_name:
                    result_dict[experiment_id] = {'name': experiment_name,
                                                  'date': experiment_date,
                                                   'progress': experiment_progress}

    return result_dict


def _experiment_exists(base_dir, experiment_id):
    return os.path.exists(os.path.join(base_dir, experiment_id))


class ExperimentsLoader(ABC):
    @abstractmethod
    def load_experiments(self) -> Dict:
        pass

    @abstractmethod
    def experiment_exists(self, experiment_id) -> bool:
        pass

    @abstractmethod
    def store_experiment(self, experiment_id: str, result, filename: str = None):
        pass

    def _delete_experiment(self, base_dir, experiment_id):
        shutil.rmtree(os.path.join(base_dir, experiment_id), ignore_errors=True)

    def _rename_experiment(self, base_dir, experiment_id, experiment_name):
        metadata_path = os.path.join(base_dir, experiment_id, "metadata.json")
        print(metadata_path)
        with open(metadata_path, 'r') as f:
            metadata = json.load(f)
        metadata['name'] = experiment_name
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=4)


class ConformationsExperimentsLoader(ExperimentsLoader):
    def __init__(self):
        self.conformations_file_name = 'conformations.pdb'

    def experiment_exists(self, experiment_id) -> bool:
        return _experiment_exists(CONFORMATIONS_EXPERIMENTS_DIR, experiment_id)

    def load_experiments(self) -> Dict:
        return _load_experiments_ids_names(CONFORMATIONS_EXPERIMENTS_DIR)

    def load_experiment(self, experiment_id) -> Dict:
        result = {}
        experiment_dir = os.path.join(CONFORMATIONS_EXPERIMENTS_DIR, experiment_id)
        loader = FileLoaderFactory().get_loader(self.conformations_file_name)
        loaded_content = loader.load(experiment_dir, self.conformations_file_name)
        result['pdb'] = loaded_content
        return result

    def store_experiment(self, experiment_id: str, result, filename: str = None):
        file_saver = FileSaverFactory().get_saver(self.conformations_file_name)
        experiment_dir = os.path.join(CONFORMATIONS_EXPERIMENTS_DIR, experiment_id)
        file_saver.save(result, experiment_dir, self.conformations_file_name)

    def save_experiment_metadata(self, experiment_id: str, experiment_name: str):
        metadata = {
            "id": experiment_id,
            "name": experiment_name,
            "date": datetime.now().isoformat()
        }

        metadata_path = os.path.join(CONFORMATIONS_EXPERIMENTS_DIR, experiment_id, "metadata.json")
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=4)

    def read_experiment_metadata(self, experiment_id: str):
        metadata_path = os.path.join(CONFORMATIONS_EXPERIMENTS_DIR, experiment_id, "metadata.json")

        # Check if metadata.json exists
        if not os.path.exists(metadata_path):
            raise FileNotFoundError(f"No metadata found for experiment_id: {experiment_id}")

        with open(metadata_path, 'r') as f:
            metadata = json.load(f)

        return metadata

    def delete_experiment(self, experiment_id):
        self._delete_experiment(CONFORMATIONS_EXPERIMENTS_DIR, experiment_id)

    def rename_experiment(self, experiment_id, experiment_name):
        self._rename_experiment(CONFORMATIONS_EXPERIMENTS_DIR, experiment_id, experiment_name)


class ProteinLabExperimentsLoader(ExperimentsLoader):
    def __init__(self):
        self.task2results_map = {
            "folding": "folding_prediction.pdb",
            "gene_ontology": "gene_ontology.json",
            "localisation": "cell_localisation.json",
            "solubility": "solubility.json"
        }
        self.fasta_saver = FastaFileSaver()

    def store_single_input(self, content, filename, experiment_id):
        inputs_dir = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id, 'inputs')
        self.fasta_saver.save(content, inputs_dir, filename)

    def store_inputs(self, fasta_files: List[FileStorage], experiment_id):
        inputs_dir = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id, 'inputs')

        for fasta_file in fasta_files:
            self.fasta_saver.save(fasta_file, inputs_dir, fasta_file.filename)

    def get_input_files(self, experiment_id):
        inputs_dir = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id, 'inputs')
        all_sequence_files = []

        # Traverse the experiment directory
        for root, dirs, files in os.walk(inputs_dir):
            for file in files:
                if file.endswith('.fasta'):
                    file_path = os.path.join(root, file)
                    all_sequence_files.append(file_path)

        return all_sequence_files
    
    def get_sequences(self, file_paths):
        fasta_file_loader = FastaFileLoader()
        all_sequence_ids = []
        all_sequences = []

        for file_path in file_paths:
            sequence_ids, sequences = fasta_file_loader.load(file_path=file_path)
            all_sequence_ids.extend(sequence_ids)
            all_sequences.extend(sequences)

        return all_sequence_ids, all_sequences

    def experiment_exists(self, experiment_id) -> bool:
        return _experiment_exists(PROTEIN_EXPERIMENTS_DIR, experiment_id)

    def load_experiments(self) -> Dict:
        return _load_experiments_ids_names(PROTEIN_EXPERIMENTS_DIR)
    
    def get_protein_ids(self, experiment_id):
        input_fastas = self.get_input_files(experiment_id)
        protein_ids, _ = self.get_sequences(input_fastas)
        return protein_ids

    def load_predictions(self, experiment_id, protein_id) -> Dict:
        result = {}
        for key, filename in self.task2results_map.items():
            output_dir = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id, protein_id)
            loader = FileLoaderFactory().get_loader(filename)
            loaded_content = loader.load(output_dir, filename)
            result[key] = loaded_content
        return result
    
    def load_experiment_progress(self, experiment_id):
        experiment_dir = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id)
        return get_progress(experiment_dir)
    
    def load_protein_progress(self, experiment_id, protein_id):
        experiment_dir = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id)
        res =  get_protein_progress(experiment_dir, protein_id)
        print("Protein_progress: ", res)
        return res
    
    def save_outputs(self, result, task, save_dir):
        filename = self.task2results_map[task]
        file_saver = FileSaverFactory().get_saver(filename)
        file_saver.save(result, save_dir, filename)

    def store_experiment(self, experiment_id: str, sequence, result, task):
        save_dir = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id, sequence)
        self.save_outputs(result, task, save_dir)

    def save_experiment_metadata(self, experiment_id: str, experiment_name: str):
        metadata = {
            "id": experiment_id,
            "name": experiment_name,
            "date": datetime.now().isoformat()
        }

        metadata_path = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id, "metadata.json")
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=4)

    def read_experiment_metadata(self, experiment_id: str):
        metadata_path = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id, "metadata.json")

        # Check if metadata.json exists
        if not os.path.exists(metadata_path):
            raise FileNotFoundError(f"No metadata found for experiment_id: {experiment_id}")

        with open(metadata_path, 'r') as f:
            metadata = json.load(f)

        return metadata

    def delete_experiment(self, experiment_id):
        self._delete_experiment(PROTEIN_EXPERIMENTS_DIR, experiment_id)

    def rename_experiment(self, experiment_id, experiment_name):
        self._rename_experiment(PROTEIN_EXPERIMENTS_DIR, experiment_id, experiment_name)


class DTILabExperimentsLoader(ExperimentsLoader):
    def __init__(self):
        self.task2results_map = {
            "dti": "skip",
        }
        self.loader = DTILoader()

    def experiment_exists(self, experiment_id) -> bool:
        return _experiment_exists(DTI_EXPERIMENTS_DIR, experiment_id)

    def load_experiments(self) -> Dict:
        return _load_experiments_ids_names(DTI_EXPERIMENTS_DIR)

    def load_result(self, experiment_id: str, protein_id: str, ligand_id: str):
        loaded_content = self.loader.get_dti_single_result(DTI_EXPERIMENTS_DIR, experiment_id, protein_id, ligand_id)
        return loaded_content
    
    def get_protein_ids(self, experiment_id: str):
        return self.loader.get_protein_ids(experiments_folder=DTI_EXPERIMENTS_DIR, experiment_id=experiment_id)
    
    def get_ligands_ids(self, experiment_id: str, protein_id: str):
        protein_folder = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id, protein_id, 'result')
        return self.loader.get_ligand_names(protein_folder)
    
    def load_experiment_progress(self, experiment_id):
        experiment_dir = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id)
        return get_progress(experiment_dir)

    def save_experiment_metadata(self, experiment_id: str, experiment_name: str):
        metadata = {
            "id": experiment_id,
            "name": experiment_name,
            "date": datetime.now().isoformat()
        }

        metadata_path = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id, "metadata.json")
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=4)

    def delete_experiment(self, experiment_id):
        self._delete_experiment(DTI_EXPERIMENTS_DIR, experiment_id)

    def rename_experiment(self, experiment_id, experiment_name):
        self._rename_experiment(DTI_EXPERIMENTS_DIR, experiment_id, experiment_name)

    def read_experiment_metadata(self, experiment_id: str):
        metadata_path = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id, "metadata.json")

        # Check if metadata.json exists
        if not os.path.exists(metadata_path):
            raise FileNotFoundError(f"No metadata found for experiment_id: {experiment_id}")

        with open(metadata_path, 'r') as f:
            metadata = json.load(f)

        return metadata

    def store_experiment(self, experiment_id: str, result, filename: str):
        file_saver = FileSaverFactory().get_saver(filename)
        experiment_dir = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id)
        file_saver.save(result, experiment_dir, filename)
