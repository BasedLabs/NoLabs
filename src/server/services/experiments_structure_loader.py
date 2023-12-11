import json
import os
import shutil
from abc import ABC, abstractmethod
from datetime import datetime
from typing import List, Dict

from src.server.services.loaders import FileLoaderFactory, DTILoader, ProteinDesignLoader
from src.server.services.savers import FileSaverFactory
from src.server.settings import PROTEIN_EXPERIMENTS_DIR, DTI_EXPERIMENTS_DIR, CONFORMATIONS_EXPERIMENTS_DIR, \
    PROTEIN_DESIGN_EXPERIMENTS_DIR


class ExperimentsLoader(ABC):
    def __init__(self, experiments_dir):
        self.experiments_dir = experiments_dir

    def load_experiments(self) -> Dict:
        experiment_ids = [d for d in os.listdir(self.experiments_dir) if
                          os.path.isdir(os.path.join(self.experiments_dir, d))]

        experimentId2name = {}

        for experiment_id in experiment_ids:
            metadata_path = os.path.join(self.experiments_dir, experiment_id, 'metadata.json')

            # Check if metadata.json exists in the directory
            if os.path.exists(metadata_path):
                with open(metadata_path, 'r') as f:
                    data = json.load(f)
                    # Assuming the JSON structure has a key 'experiment_name' with the name of the experiment.
                    # Adjust if the structure is different.
                    experiment_name = data.get('name', None)
                    if experiment_name:
                        experimentId2name[experiment_id] = experiment_name

        return experimentId2name

    def experiment_exists(self, experiment_id) -> bool:
        return os.path.exists(os.path.join(self.experiments_dir, experiment_id))

    @abstractmethod
    def store_experiment(self, experiment_id: str, result, filename: str = None):
        pass

    @abstractmethod
    def load_experiment(self, experiment_id) -> Dict:
        pass

    def delete_experiment(self, experiment_id):
        shutil.rmtree(os.path.join(self.experiments_dir, experiment_id), ignore_errors=True)

    def rename_experiment(self, experiment_id, experiment_name):
        metadata_path = os.path.join(self.experiments_dir, experiment_id, "metadata.json")
        print(metadata_path)
        with open(metadata_path, 'r') as f:
            metadata = json.load(f)
        metadata['name'] = experiment_name
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=4)

    def save_experiment_metadata(self, experiment_id: str, experiment_name: str):
        metadata = {
            "id": experiment_id,
            "name": experiment_name,
            "date": datetime.now().isoformat()
        }

        metadata_path = os.path.join(self.experiments_dir, experiment_id, "metadata.json")
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=4)

    def read_experiment_metadata(self, experiment_id: str):
        metadata_path = os.path.join(self.experiments_dir, experiment_id, "metadata.json")

        # Check if metadata.json exists
        if not os.path.exists(metadata_path):
            raise FileNotFoundError(f"No metadata found for experiment_id: {experiment_id}")

        with open(metadata_path, 'r') as f:
            metadata = json.load(f)

        return metadata


class ConformationsExperimentsLoader(ExperimentsLoader):
    def __init__(self):
        super().__init__(CONFORMATIONS_EXPERIMENTS_DIR)
        self.conformations_file_name = 'conformations.pdb'

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


class ProteinLabExperimentsLoader(ExperimentsLoader):
    def __init__(self):
        super().__init__(PROTEIN_EXPERIMENTS_DIR)
        self.task2results_map = {
            "folding": "folding_prediction.pdb",
            "gene_ontology": "gene_ontology.json",
            "localisation": "cell_localisation.json",
            "solubility": "solubility.json"
        }

    def load_experiment(self, experiment_id) -> Dict:
        result = {}
        for key, filename in self.task2results_map.items():
            experiment_dir = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id)
            loader = FileLoaderFactory().get_loader(filename)
            loaded_content = loader.load(experiment_dir, filename)
            result[key] = loaded_content
        return result

    def store_experiment(self, experiment_id: str, result, filename: str = None):
        file_saver = FileSaverFactory().get_saver(filename)
        experiment_dir = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id)
        file_saver.save(result, experiment_dir, filename)

    def read_experiment_metadata(self, experiment_id: str):
        metadata_path = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id, "metadata.json")

        # Check if metadata.json exists
        if not os.path.exists(metadata_path):
            raise FileNotFoundError(f"No metadata found for experiment_id: {experiment_id}")

        with open(metadata_path, 'r') as f:
            metadata = json.load(f)

        return metadata


class DTILabExperimentsLoader(ExperimentsLoader):
    def __init__(self):
        super().__init__(DTI_EXPERIMENTS_DIR)
        self.task2results_map = {
            "dti": "skip",
        }

    def load_experiment(self, experiment_id: str):
        loader = DTILoader()
        loaded_content = loader.get_dti_results(DTI_EXPERIMENTS_DIR, experiment_id)
        return loaded_content

    def read_experiment_metadata(self, experiment_id: str):
        metadata_path = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id, "metadata.json")

        # Check if metadata.json exists
        if not os.path.exists(metadata_path):
            raise FileNotFoundError(f"No metadata found for experiment_id: {experiment_id}")

        with open(metadata_path, 'r') as f:
            metadata = json.load(f)

        return metadata

    def store_experiment(self, experiment_id: str, result, filename: str = None):
        file_saver = FileSaverFactory().get_saver(filename)
        experiment_dir = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id)
        file_saver.save(result, experiment_dir, filename)


class ProteinDesignExperimentsLoader(ExperimentsLoader):
    def __init__(self):
        super().__init__(PROTEIN_DESIGN_EXPERIMENTS_DIR)
        self.result_folder = 'generated'

    def load_experiment(self, experiment_id) -> Dict:
        loader = ProteinDesignLoader()
        pdb_contents = loader.get_protein_design_results(
            PROTEIN_DESIGN_EXPERIMENTS_DIR,
            experiment_id,
            self.result_folder)
        return pdb_contents

    def store_experiment(self, experiment_id: str, result, filename: str = None):
        file_saver = FileSaverFactory().get_saver(filename)
        experiment_dir = os.path.join(self.experiments_dir, experiment_id)
        file_saver.save(result, experiment_dir, filename)
