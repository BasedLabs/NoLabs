import os
from abc import abstractmethod, ABC
from typing import List

from werkzeug.datastructures import FileStorage

from src.server.services.conformations.conformations_pipeline import permute_simulation
from src.server.services.experiments_structure_loader import ProteinLabExperimentsLoader, DTILabExperimentsLoader
from src.server.services.inference_service import create_pipeline, create_model, get_models_from_config
from src.server.services.savers import FastaFileSaver, SDFFileSaver, PDBFileSaver
from src.server.settings import EXPERIMENTS_DIR, PROTEIN_EXPERIMENTS_DIR, DTI_EXPERIMENTS_DIR, \
    CONFORMATIONS_EXPERIMENTS_DIR, PROTEIN_DESIGN_EXPERIMENTS_DIR
from src.server.services.mixins import UUIDGenerator


# Helper function to ensure the base directory exists
def ensure_base_directory():
    if not os.path.exists(EXPERIMENTS_DIR):
        os.makedirs(EXPERIMENTS_DIR)
    if not os.path.exists(PROTEIN_EXPERIMENTS_DIR):
        os.makedirs(PROTEIN_EXPERIMENTS_DIR)
    if not os.path.exists(DTI_EXPERIMENTS_DIR):
        os.makedirs(DTI_EXPERIMENTS_DIR)
    if not os.path.exists(CONFORMATIONS_EXPERIMENTS_DIR):
        os.makedirs(CONFORMATIONS_EXPERIMENTS_DIR)


ensure_base_directory()


class BaseExperiment(ABC, UUIDGenerator):
    def __init__(self):
        self.pipeline = create_pipeline(use_gpu, is_test)


    @abstractmethod
    def run(self, *args, **kwargs):
        pass


class ProteinPropertyPrediction(BaseExperiment):

    def __init__(self, use_gpu=False, is_test=True):
        super().__init__()
        self.loader = ProteinLabExperimentsLoader()
        models_metadata = get_models_from_config(is_test)
        for model_metadata in models_metadata:
            if model_metadata["task"] in self.loader.task2results_map.keys():
                model = create_model(model_metadata, use_gpu)
                self.pipeline.add_model(model)

    def run(self, sequence: str, experiment_id: str):
        if not experiment_id:
            experiment_id = self.gen_uuid()
        for task, result_file in self.loader.task2results_map.items():
            self.run_by_task(task, sequence, experiment_id, result_file)
        return experiment_id

    def run_by_task(self, task: str, sequence: str, experiment_id: str, result_file: str):
        model = self.pipeline.get_model_by_task(task)
        result = model.predict(sequence)
        self.loader.store_experiment(experiment_id, result, result_file)


class DrugDiscovery(BaseExperiment):

    def __init__(self, use_gpu, is_test=True):
        # "skip" file fornat is used for skipping since the data saving procedure is written somewhere else
        # Here it's written in dti model predict method
        super().__init__()
        self.loader = DTILabExperimentsLoader()
        models_metadata = get_models_from_config(is_test)
        for model_metadata in models_metadata:
            if model_metadata["task"] in self.loader.task2results_map.keys():
                model = create_model(model_metadata, use_gpu)
                self.pipeline.add_model(model)

    def _store_inputs(self, experiments_dir: str, ligand_files: List[FileStorage], protein_files: List[FileStorage]):
        sdf_saver = SDFFileSaver()
        pdb_saver = PDBFileSaver()
        fasta_saver = FastaFileSaver()

        ligand_file_paths = []
        pdb_file_paths = []

        for ligand_file in ligand_files:
            file_path = sdf_saver.save(ligand_file, experiments_dir, ligand_file.filename)
            ligand_file_paths.append(file_path)

        for protein_file in protein_files:
            if protein_file.filename.endswith(".pdb"):
                file_path = pdb_saver.save(protein_file, experiments_dir, protein_file.filename)
                pdb_saver.pdb_to_fasta(protein_file, experiments_dir, protein_file.filename)
                pdb_file_paths.append(file_path)
            elif protein_file.filename.endswith(".fasta"):
                file_path = fasta_saver.save(protein_file, experiments_dir, protein_file.filename)
                pdb_file_paths.append(file_path)

        return (ligand_file_paths, pdb_file_paths)

    def run(self, ligand_files: List[FileStorage], protein_files: List[FileStorage], experiment_id: str):
        if not experiment_id:
            experiment_id = self.gen_uuid()
        model = self.pipeline.get_model_by_task("dti")
        experiment_dir = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id)
        ligand_file_paths, pdb_file_paths = self._store_inputs(experiment_dir, ligand_files, protein_files)
        model.predict(ligand_file_paths, pdb_file_paths, experiment_dir)
        return experiment_id


class ProteinDesign(BaseExperiment):
    def __init__(self):
        super().__init__()

    def run(self,
            experiment_id: str,
            pdb_file_name: str = None,
            pdb_content: str = None,
            contig: str = '50',
            symmetry: str = None,
            timesteps: int = 50,
            hotspots: str = ''):
        assert pdb_file_name and pdb_content or not (pdb_file_name and pdb_content)
        if not experiment_id:
            experiment_id = self.gen_uuid()
        model = self.pipeline.get_model_by_task("protein_design")
        experiment_dir = os.path.join(PROTEIN_DESIGN_EXPERIMENTS_DIR, experiment_id)
        if pdb_file_name:
            with open(os.path.join(experiment_dir, pdb_file_name), 'w') as pdb_file:
                pdb_file.write(pdb_content)
        result = model.predict(pdb_content, contig, symmetry, timesteps, hotspots)
        return experiment_id