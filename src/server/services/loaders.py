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

    def get_dti_results(pipeline, ligand_files: str):
        model = pipeline.get_model_by_task("dti")
        experiment_folder = model.experiment_folder
        result_folder = model.result_folder
        protein_file = model.protein_file
        protein_name = model.protein_name

        pdb_content = ""

        # Open and read the PDB file
        with open(f'{experiment_folder}/{protein_file}', 'r') as pdb_file:
            for line in pdb_file:
                pdb_content += line

        ligands_sdf_contents = []
        affinity_list = []

        ligand_names = [os.path.splitext(file.filename)[0] for file in ligand_files]

        for ligand_name in ligand_names:

            ligand_file = f'{result_folder}/{ligand_name}_tankbind.sdf'

            info_df = pd.read_csv(f"{result_folder}/{ligand_name}_info_with_predicted_affinity.csv")
            chosen = info_df.loc[info_df.groupby(['protein_name', 'compound_name'],sort=False)['affinity'].agg('idxmax')].reset_index()
            affinity_list.append(chosen['affinity'].item())

            sdf_supplier = SDMolSupplier(ligand_file)
            # Initialize an empty string to store the SDF contents
            sdf_contents = ""
            # Iterate through the molecules in the SDF file and append their representations to the string
            for mol in sdf_supplier:
                if mol is not None:
                    # Convert the molecule to an SDF block and append it to the string
                    sdf_contents += Chem.MolToMolBlock(mol) + "\n"
            ligands_sdf_contents.append(sdf_contents)

        
        print("AFFINITY_LIST: ", affinity_list)

        return pdb_content, protein_name, ligands_sdf_contents, ligand_names, affinity_list
