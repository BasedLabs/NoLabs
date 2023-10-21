import os
import json
import csv

import pandas as pd
from rdkit import Chem
from rdkit.Chem import SDMolSupplier


class FileLoader:
    def load(self, folder, filename):
        raise NotImplementedError

class PDBFileLoader(FileLoader):
    def load(self, folder, filename):
        with open(os.path.join(folder, filename), 'r') as f:
            return f.read()

class CSVFileLoader(FileLoader):
    def load(self, folder, filename):
        content = []
        with open(os.path.join(folder, filename), 'r', newline='') as f:
            reader = csv.reader(f)
            for row in reader:
                content.append(row)
        return content

class SDFFileLoader(FileLoader):
    def load(self, folder, filename):
        with open(os.path.join(folder, filename), 'r') as f:
            return f.read()

class JSONFileLoader(FileLoader):
    def load(self, folder, filename):
        with open(os.path.join(folder, filename), 'r') as f:
            return json.load(f)

class CSVListLoader(FileLoader):
    def load(self, folder, filename):
        content = []
        with open(os.path.join(folder, filename), 'r', newline='') as f:
            reader = csv.reader(f)
            for row in reader:
                row[1] = float(row[1])  # Convert second element to float
                content.append(row)
        return content

class FileLoaderFactory:
    @staticmethod
    def get_loader(filename):
        file_extension = os.path.splitext(filename)[1].lower()
        if file_extension == '.pdb':
            return PDBFileLoader()
        elif file_extension == '.csv':
            return CSVFileLoader()
        elif file_extension == '.sdf':
            return SDFFileLoader()
        elif file_extension == '.json':
            return JSONFileLoader()
        elif file_extension == '.list.csv':
            return CSVListLoader()
        else:
            raise ValueError(f"Unsupported file extension: {file_extension}")

class DTILoader:
    def __init__(self):
        pass

    def get_dti_results(self, experiment_folder: str):
        protein_names = [d for d in os.listdir(experiment_folder) \
        if os.path.isdir(os.path.join(experiment_folder, d))]

        results = []
        for protein_name in protein_names:
            result_folder = f'{experiment_folder}{protein_name}/result'
            protein_name = protein_name

            pdb_content = ""

            # Open and read the PDB file
            with open(f'{experiment_folder}/{protein_name}/{protein_name}.pdb', 'r') as pdb_file:
                for line in pdb_file:
                    pdb_content += line

            ligand_names = self.get_ligand_names(f'{experiment_folder}/{protein_name}')

            for ligand_name in ligand_names:

                ligand_file = f'{result_folder}/{ligand_name}_tankbind.sdf'

                info_df = pd.read_csv(f"{result_folder}/{ligand_name}_info_with_predicted_affinity.csv")
                chosen = info_df.loc[info_df.groupby(['protein_name', 'compound_name'],sort=False)['affinity'].agg('idxmax')].reset_index()

                sdf_supplier = SDMolSupplier(ligand_file)
                # Initialize an empty string to store the SDF contents
                sdf_contents = ""
                # Iterate through the molecules in the SDF file and append their representations to the string
                for mol in sdf_supplier:
                    if mol is not None:
                        # Convert the molecule to an SDF block and append it to the string
                        sdf_contents += Chem.MolToMolBlock(mol) + "\n"
                results.append(
                    {
                    'pdb': pdb_content, 
                    'proteinName': protein_name, 
                    'sdf': sdf_contents, 
                    'ligandName': ligand_name, 
                    'affinity': chosen['affinity'].item()
                    }
                )

        return results

    def get_ligand_names(self, ligands_path):
        # List all directories inside the ligands folder
        all_dirs = [d for d in os.listdir(ligands_path) if os.path.isdir(os.path.join(ligands_path, d))]

        # Extract ligand names by removing the '_dataset' suffix
        ligand_names = [dir_name.rsplit('_dataset', 1)[0] for dir_name in all_dirs if dir_name.endswith('_dataset')]

        return ligand_names


def load_experiment_names(experiments_dir):
    # List all directories inside the dti experiments folder
    experiment_ids = [d for d in os.listdir(experiments_dir) if os.path.isdir(os.path.join(experiments_dir, d))]

    experimentId2name = {}

    for experiment_id in experiment_ids:
        metadata_path = os.path.join(experiments_dir, experiment_id, 'metadata.json')

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
