import os
import shutil
from abc import abstractmethod, ABC
from typing import List
import numpy as np

from werkzeug.datastructures import FileStorage

from nolabs.server.services.experiments_structure_loader import ProteinLabExperimentsLoader, DTILabExperimentsLoader
from nolabs.server.services.inference_service import create_pipeline, create_model, get_models_from_config
from nolabs.server.services.savers import FastaFileSaver, SDFFileSaver, PDBFileSaver
from nolabs.server.settings import EXPERIMENTS_DIR, PROTEIN_EXPERIMENTS_DIR, DTI_EXPERIMENTS_DIR, \
    CONFORMATIONS_EXPERIMENTS_DIR, PROTEIN_DESIGN_EXPERIMENTS_DIR
from nolabs.server.services.mixins import UUIDGenerator

from nolabs.ai.model import PocketPredictor


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
        self.loader = ProteinLabExperimentsLoader()
        self.pipeline = create_pipeline(use_gpu, is_test, target_tasks=["localisation",
                                                                         "folding",
                                                                         "solubility",
                                                                         "gene_ontology"])

    def run(self, sequence: str, protein_files: List[FileStorage], experiment_id: str):

        if not experiment_id:
            experiment_id = self.gen_uuid()

        if sequence:
            fasta_contents = f">{sequence}\n{sequence}"
            filename = f"{sequence}.fasta"
            self.loader.store_single_input(fasta_contents, sequence, experiment_id=experiment_id)
            for task, _ in self.loader.task2results_map.items():
                self.run_single_expeeriment(task, sequence, experiment_id)

        if protein_files:
            self.loader.store_inputs(fasta_files=protein_files, experiment_id=experiment_id)
            input_fastas = self.loader.get_input_files(experiment_id=experiment_id)
            experiment_dir = os.path.join(PROTEIN_EXPERIMENTS_DIR, experiment_id)
            self.pipeline.run(input_fastas, loader=self.loader, experiment_dir=experiment_dir)

        return experiment_id

    def run_single_experiment(self, task: str, sequence: List[str], experiment_id: str):
        model = self.pipeline.get_model_by_task(task)
        result = model.predict(sequence)
        self.loader.store_experiment(experiment_id, sequence, result, task)


class DrugDiscovery(BaseExperiment):

    def __init__(self, use_gpu, is_test=True):
        # "skip" file fornat is used for skipping since the data saving procedure is written somewhere else
        # Here it's written in dti model predict method
        self.loader = DTILabExperimentsLoader()
        self.pipeline = create_pipeline(use_gpu, is_test, target_tasks=["dti"])

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
                pdb_saver.pdb_to_fasta(protein_file.stream, experiments_dir, protein_file.filename)
                pdb_file_paths.append(file_path)
            elif protein_file.filename.endswith(".fasta"):
                file_path = fasta_saver.save(protein_file, experiments_dir, protein_file.filename)
                pdb_file_paths.append(file_path)

        return (ligand_file_paths, pdb_file_paths)

    def run(self, experiment_id: str):
        if not experiment_id:
            experiment_id = self.gen_uuid()
        model = self.pipeline.get_model_by_task("dti")
        experiment_dir = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id)

        binding_pockets = []
        protein_file_paths = []
        protein_ids = []
        ligand_file_paths = [] 

        targets_dir = os.path.join(experiment_dir, 'targets')
           # List all items in the targets_dir
        for item in os.listdir(targets_dir):
            subdir = os.path.join(targets_dir, item)
            
            # Check if the item is a directory
            if os.path.isdir(subdir):
                target_id = os.path.basename(subdir)
                protein_ids.append(target_id)
                pocket_dir = os.path.join(subdir, 'pocket')
                pocket_file = os.path.join(pocket_dir, 'pocket.npy')

                if os.path.isdir(pocket_dir) and os.path.isfile(pocket_file):
                    # Load the numpy array if pocket.npy exists
                    binding_pocket = np.load(pocket_file)
                    binding_pockets.append(binding_pocket)
                    for file in os.listdir(subdir):
                        if file.endswith('.fasta'):
                            pdb_file_path = os.path.join(subdir, file)
                            protein_file_paths.append(pdb_file_path)
                else:
                    # Check for .pdb files in the subdir
                    for file in os.listdir(subdir):
                        if file.endswith('.fasta'):
                            pdb_file_path = os.path.join(subdir, file)
                            protein_file_paths.append(pdb_file_path)
                            self.loader.predict_3d_structure(experiment_id, target_id)
                            # Assuming predict_pocket is a method that predicts the pocket and returns it
                            binding_pocket = self.predict_pocket(experiment_id=experiment_id, protein_id=target_id)
                            binding_pockets.append(np.array(binding_pocket))

        ligands_dir = os.path.join(experiment_dir, 'ligands')
        for subdir, dirs, files in os.walk(ligands_dir):
            for file in files:
                if file.endswith('.sdf'):
                    ligand_file_paths.append(os.path.join(subdir, file))

                            
        print("BINDING POCKETS: ", binding_pockets)
        print("Ligand file paths: ", ligand_file_paths)
        print("Protein file paths: ", protein_file_paths)
        print("Protein Ids:", protein_ids)

        model.predict(ligand_file_paths, protein_file_paths, protein_ids, binding_pockets, experiment_dir)

        return experiment_id
    

    def predict_pocket(self, experiment_id: str, protein_id: str):

        pocket_predictor = PocketPredictor('pocket_predictor', 'pocket_prediction')
        pocket_predictor.load_model()

        protein_dir = os.path.join(DTI_EXPERIMENTS_DIR, experiment_id, 'targets', protein_id)

        pockets_dir = os.path.join(protein_dir, 'pocket')
        if not os.path.exists(pockets_dir):
            os.mkdir(pockets_dir)

        protein_file = None

        for filename in os.listdir(protein_dir):
            if filename.endswith(".pdb"):
                destination_protein_file = os.path.join(pockets_dir, filename)
                shutil.copyfile(os.path.join(protein_dir, filename), destination_protein_file)
                protein_file = os.path.join(pockets_dir, filename)

        return pocket_predictor.predict(protein_file, pockets_dir)



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