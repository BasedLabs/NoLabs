import json
import os
import glob
import shutil
from abc import ABC, abstractmethod
from dataclasses import dataclass
from datetime import datetime
from typing import List, Dict
from pathlib import Path
from werkzeug.datastructures import FileStorage
from pathlib import Path
from typing import List, Dict, Tuple
import numpy as np

from io import StringIO

from src.server.services.loaders import FileLoaderFactory, DTILoader, FastaFileLoader,PDBFileLoader, SDFFileLoader
from src.server.services.savers import FileSaverFactory, FastaFileSaver, SDFFileSaver, PDBFileSaver
from src.server.settings import PROTEIN_EXPERIMENTS_DIR, DTI_EXPERIMENTS_DIR, CONFORMATIONS_EXPERIMENTS_DIR
from src.server.services.progress import get_progress
from src.server.services.mixins import UUIDGenerator
from src.server.services.loaders import FileLoaderFactory, DTILoader
from src.server.services.savers import FileSaverFactory
from src.server.settings import PROTEIN_EXPERIMENTS_DIR, DTI_EXPERIMENTS_DIR, CONFORMATIONS_EXPERIMENTS_DIR, \
    PROTEIN_DESIGN_EXPERIMENTS_DIR
from test.ai.mock_model import APIFolding


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
        target_dir =  os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id, protein_id)
        res =  get_progress(target_dir)
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


class ProteinDesignExperimentsLoader(ExperimentsLoader):
    def experiment_exists(self, experiment_id) -> bool:
        return _experiment_exists(PROTEIN_DESIGN_EXPERIMENTS_DIR, experiment_id)

    def load_experiments(self) -> Dict:
        return _load_experiments_ids_names(PROTEIN_DESIGN_EXPERIMENTS_DIR)

    def load_experiment(self, experiment_id) -> List[str]:
        experiment_dir = os.path.join(PROTEIN_DESIGN_EXPERIMENTS_DIR, experiment_id)
        pdbs = []
        for file in glob.glob(os.path.join(experiment_dir, '*_designed*.pdb')):
            with open(file, 'r') as f:
                pdbs.append(f.read())
        return pdbs

    def store_experiment(self, experiment_id: str, result: List[str], filename: str = None):
        for i, pdb in enumerate(result):
            fn_without_extension, extension = os.path.splitext(filename)
            fn_without_extension += '_designed_' + str(i)
            fn = fn_without_extension + extension
            file_saver = FileSaverFactory().get_saver(fn)
            experiment_dir = os.path.join(PROTEIN_DESIGN_EXPERIMENTS_DIR, experiment_id)
            file_saver.save(pdb, experiment_dir, fn)

    def save_experiment_metadata(self, experiment_id: str, experiment_name: str):
        metadata = {
            "id": experiment_id,
            "name": experiment_name,
            "date": datetime.now().isoformat()
        }

        metadata_path = os.path.join(PROTEIN_DESIGN_EXPERIMENTS_DIR, experiment_id, "metadata.json")
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=4)

    def read_experiment_metadata(self, experiment_id: str):
        metadata_path = os.path.join(PROTEIN_DESIGN_EXPERIMENTS_DIR, experiment_id, "metadata.json")

        # Check if metadata.json exists
        if not os.path.exists(metadata_path):
            raise FileNotFoundError(f"No metadata found for experiment_id: {experiment_id}")

        with open(metadata_path, 'r') as f:
            metadata = json.load(f)

        return metadata

    def delete_experiment(self, experiment_id):
        self._delete_experiment(PROTEIN_DESIGN_EXPERIMENTS_DIR, experiment_id)

    def rename_experiment(self, experiment_id, experiment_name):
        self._rename_experiment(PROTEIN_DESIGN_EXPERIMENTS_DIR, experiment_id, experiment_name)


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

    def check_result_available(self, experiment_id: str, protein_id: str, ligand_id: str):
        return self.loader.check_result_available(DTI_EXPERIMENTS_DIR, experiment_id, protein_id, ligand_id)

    def get_protein_ids(self, experiment_id: str):
        return self.loader.get_protein_ids(experiments_folder=DTI_EXPERIMENTS_DIR, experiment_id=experiment_id)

    def get_ligands_ids(self, experiment_id: str, protein_id: str):
        protein_folder = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id, 'results', protein_id)
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

    def save_target_metadata(self, experiment_id: str, target_id: str, protein_name: str, manual_pocket = False):
        metadata = {
            "id": target_id,
            "name": protein_name,
            "date": datetime.now().isoformat(),
            "manual_pocket": manual_pocket
        }

        metadata_path = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id, 'targets', target_id, "metadata.json")
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=4)

    def save_ligand_metadata(self, experiment_id: str, ligand_id: str, ligand_name: str):
        metadata = {
            "id": ligand_id,
            "name": ligand_name,
            "date": datetime.now().isoformat(),
        }

        metadata_path = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id, 'ligands', ligand_id, "metadata.json")
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=4)

    def update_target_metadata(self, experiment_id: str, target_id: str, key: str, value):
        metadata = self.load_target_metadata(experiment_id, target_id)
        metadata[key] = value

        metadata_path = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id, 'targets', target_id, "metadata.json")
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=4)


    def load_target_metadata(self, experiment_id: str, target_id: str,):
        metadata_path = os.path.join(DTI_EXPERIMENTS_DIR,
                                      experiment_id,
                                        'targets',
                                          target_id,
                                            "metadata.json")

        # Check if metadata.json exists
        if not os.path.exists(metadata_path):
            return {}

        with open(metadata_path, 'r') as f:
            metadata = json.load(f)

        return metadata
    

    def load_ligand_metadata(self, experiment_id: str, ligand_id: str,):
        metadata_path = os.path.join(DTI_EXPERIMENTS_DIR,
                                      experiment_id,
                                        'ligands',
                                          ligand_id,
                                            "metadata.json")

        # Check if metadata.json exists
        if not os.path.exists(metadata_path):
            return {}

        with open(metadata_path, 'r') as f:
            metadata = json.load(f)

        return metadata

    def store_target(self, experiment_id: str, protein_file: List[FileStorage]):
        experiments_dir = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id)
        targets_dir = os.path.join(experiments_dir, 'targets')

        if not os.path.exists(targets_dir):
            os.mkdir(targets_dir)

        pdb_saver = PDBFileSaver()
        fasta_saver = FastaFileSaver()

        id_generator = UUIDGenerator()

        protein_id = id_generator.gen_uuid()['id']
        protein_name = os.path.splitext(protein_file.filename)[-2]
        protein_dir = os.path.join(targets_dir, protein_id)
        if not os.path.exists(protein_dir):
            os.mkdir(protein_dir)
        self.save_target_metadata(experiment_id, protein_id, protein_name)

        if protein_file.filename.endswith(".pdb"):
            print("Saving pdb file")
            pdb_content_str = protein_file.read().decode('utf-8')
            pdb_content = StringIO(pdb_content_str)
            pdb_saver.save(pdb_content_str, protein_dir, protein_file.filename)
            pdb_content.seek(0)
            pdb_saver.pdb_to_fasta(pdb_content, protein_dir, protein_name)

        elif protein_file.filename.endswith(".fasta"):
            fasta_saver.save(protein_file, protein_dir, protein_file.filename)

    def load_targets(self, experiment_id: str):
        targets_dir = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id, 'targets')

        pdb_loader = PDBFileLoader()
        fasta_loader = FastaFileLoader()

        targets = {}

        # Iterate over each subdirectory in the targets directory
        if not os.path.exists(targets_dir):
            return targets

        for target_id in os.listdir(targets_dir):
            target_path = os.path.join(targets_dir, target_id)

            # Skip if not a directory
            if not os.path.isdir(target_path):
                continue

            # Initialize the target entry
            targets[target_id] = {
                'metadata': self.load_target_metadata(experiment_id, target_id),
                'pdb': None,
                'fasta': None
            }

            # Check for PDB and FASTA files
            for file in os.listdir(target_path):
                if file.endswith('.pdb'):
                    file = pdb_loader.load(target_path, os.path.basename(file))
                    targets[target_id]['pdb'] = file  # Storing file name for simplicity
                elif file.endswith('.fasta'):
                    file = fasta_loader.load(os.path.join(targets_dir, target_id, file))[1][0]
                    targets[target_id]['fasta'] = file  # Storing file name for simplicity

        return targets

    def predict_3d_structure(self, experiment_id, protein_id):
        protein_dir = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id, 'targets', protein_id)

        # Find a .fasta file in the directory
        fasta_file = None
        for file in os.listdir(protein_dir):
            if file.endswith(".fasta"):
                fasta_file = file
                break

        if fasta_file is None:
            print("No .fasta file found in the directory.")
            return

        # Initialize the APIFolding model
        folding_model = APIFolding(model_name="esm_model", gpu=False)

        # Read the fasta file
        fasta_path = os.path.join(protein_dir, fasta_file)
        with open(fasta_path, 'r') as file:
            sequence = file.read()

        # Predict the pdb structure
        pdb_structure = folding_model.predict(sequence)

        # Write the prediction to <fasta name>.pdb
        pdb_file = os.path.splitext(fasta_file)[0] + '.pdb'
        pdb_path = os.path.join(protein_dir, pdb_file)
        with open(pdb_path, 'w') as file:
            file.write(pdb_structure)

        return pdb_structure
    
    def set_binding_pocket(self, experiment_id, protein_id, pocket_ids_array):
        protein_dir = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id, 'targets', protein_id)

        pockets_dir = os.path.join(protein_dir, 'pocket')
        if not os.path.exists(pockets_dir):
            os.mkdir(pockets_dir)

        np.save(os.path.join(pockets_dir, "pocket.npy"), pocket_ids_array)
        self.update_target_metadata(experiment_id, protein_id, "manual_pocket", True)


    def load_binding_pocket(self, experiment_id, protein_id):
        protein_dir = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id, 'targets', protein_id)

        pockets_dir = os.path.join(protein_dir, 'pocket')
        if not os.path.exists(pockets_dir):
            return []
        
        binding_pockets_ids = np.load(os.path.join(pockets_dir, "pocket.npy"))
        binding_pockets_ids = binding_pockets_ids.tolist()

        return binding_pockets_ids
        
    def store_ligand(self, experiment_id: str, ligand_file: List[FileStorage]):
        experiments_dir = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id)
        ligands_dir = os.path.join(experiments_dir, 'ligands')

        if not os.path.exists(ligands_dir):
            os.mkdir(ligands_dir)

        sdf_saver = SDFFileSaver()

        id_generator = UUIDGenerator()

        ligand_id = id_generator.gen_uuid()['id']
        ligand_name = os.path.splitext(ligand_file.filename)[-2]
        ligand_dir = os.path.join(ligands_dir, ligand_id)
        if not os.path.exists(ligand_dir):
            os.mkdir(ligand_dir)
        self.save_ligand_metadata(experiment_id, ligand_id, ligand_name)

        sdf_saver.save(ligand_file, ligand_dir, ligand_file.filename)
        

    def load_ligands(self, experiment_id: str):
        ligands_dir = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id, 'ligands')
        sdf_loader = SDFFileLoader()
        ligands = {}

        if not os.path.exists(ligands_dir):
            return ligands

        for ligand_id in os.listdir(ligands_dir):
            ligand_path = os.path.join(ligands_dir, ligand_id)

            if not os.path.isdir(ligand_path):
                continue

            ligands[ligand_id] = {
                'metadata': self.load_ligand_metadata(experiment_id, ligand_id),
                'sdf': None
            }

            # Check for PDB and FASTA files
            for file in os.listdir(ligand_path):
                if file.endswith('.sdf'):
                    file = sdf_loader.load(ligand_path, os.path.basename(file))
                    ligands[ligand_id]['sdf'] = file  # Storing file name for simplicity

        return ligands
    
    def get_pockets_prediction_progress(experiment_id: str):
        experiment_dir = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id)

        targets_dir = os.path.join(experiment_dir, 'targets')
        targets_with_pockets = []
        targets_without_pockets = []

        # List all items in the targets_dir
        for item in os.listdir(targets_dir):
            subdir = os.path.join(targets_dir, item)
            
            # Check if the item is a directory
            if os.path.isdir(subdir):
                target_id = os.path.basename(subdir)
                pocket_dir = os.path.join(subdir, 'pocket')
                pocket_file = os.path.join(pocket_dir, 'pocket.npy')

                if os.path.isdir(pocket_dir) and os.path.isfile(pocket_file):
                    targets_with_pockets.append(target_id)
                else:
                    targets_without_pockets.append(target_id)

        total_targets = len(targets_with_pockets) + len(targets_without_pockets)
        percentage_with_pockets = (len(targets_with_pockets) / total_targets * 100) if total_targets > 0 else 0

        return {
            'targets_with_pockets': targets_with_pockets,
            'targets_without_pockets': targets_without_pockets,
            'percentage': percentage_with_pockets
        }



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
