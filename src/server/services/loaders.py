import glob
import os
import json
import csv

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import SDMolSupplier
from pathlib import Path


class FileLoader:
    def load(self, folder, filename):
        raise NotImplementedError

class PDBFileLoader(FileLoader):
    def load(self, folder, filename):
        with open(os.path.join(folder, filename), 'r') as f:
            return f.read()

class FastaFileLoader(FileLoader):
    def load(self, file_path):
        sequence_ids = []
        sequences = []

        with open(file_path, 'r') as file:
            current_sequence_id = None
            sequence_data = []

            for line in file:
                if line.startswith('>'):
                    if current_sequence_id:
                        sequences.append(''.join(sequence_data))
                        sequence_data = []

                    current_sequence_id = line[1:].strip()
                    sequence_ids.append(current_sequence_id)
                else:
                    sequence_data.append(line.strip())

            if current_sequence_id:
                sequences.append(''.join(sequence_data))

        return sequence_ids, sequences

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
        elif file_extension == '.txt':
            return JSONFileLoader()
        elif file_extension == '.list.csv':
            return CSVListLoader()
        elif file_extension == '.fasta':
            return FastaFileLoader()
        else:
            raise ValueError(f"Unsupported file extension: {file_extension}")


class ProteinDesignLoader:
    def get_protein_design_results(self, experiments_folder: str, experiment_id: str, result_folder: str):
        import glob
        experiment_folder = os.path.join(experiments_folder, experiment_id, result_folder)
        generated_proteins_filenames = glob.glob(os.path.join(experiment_folder, '*.pdb'))
        pdb_contents = []
        for filename in generated_proteins_filenames:
            with open(filename, 'r') as f:
                pdb_contents.append(f.read())
        return pdb_contents

class DTILoader:
    def __init__(self):
        pass

    def get_dti_single_result(self,
                              experiment_folder: str,
                              experiment_id: str,
                              protein_id: str,
                              ligand_id: str):

        protein_folder = os.path.join(experiment_folder, experiment_id, 'results', protein_id)
        result_folder = os.path.join(protein_folder, ligand_id, 'result')

        pdb_content = ""
        # Open and read the PDB file
        for filename in os.listdir(result_folder):
            if filename.endswith('_pred_protein.pdb') and os.path.isfile(os.path.join(result_folder, filename)):
                pred_pdb_file_path = os.path.join(result_folder, filename)
                with open(pred_pdb_file_path, 'r') as pdb_file:
                    for line in pdb_file:
                        pdb_content += line

        ligand_file = f'{result_folder}/{ligand_id}_pred_ligand.sdf'

        plddt_df = pd.read_csv(f"{result_folder}/{ligand_id}_ligand_plddt.csv", header=None)
        ligand_plddt = np.round(plddt_df[0].mean(),1)

        sdf_supplier = SDMolSupplier(ligand_file)
        # Initialize an empty string to store the SDF contents
        sdf_contents = ""
        # Iterate through the molecules in the SDF file and append their representations to the string
        for mol in sdf_supplier:
            if mol is not None:
                # Convert the molecule to an SDF block and append it to the string
                sdf_contents += Chem.MolToMolBlock(mol) + "\n"
        return {
            'pdb': pdb_content,
            'proteinName': protein_id,
            'sdf': sdf_contents,
            'ligandName': ligand_id,
            'affinity': ligand_plddt
            }




    def get_dti_results(self, experiments_folder: str, experiment_id: str):
        experiment_folder = os.path.join(experiments_folder, experiment_id)
        protein_names = [d for d in os.listdir(experiment_folder) \
        if os.path.isdir(os.path.join(experiment_folder, d))]

        results = []
        for protein_name in protein_names:
            protein_folder = os.path.join(experiment_folder, protein_name)
            result_folder = os.path.join(protein_folder, 'result')
            protein_name = protein_name

            pdb_content = ""
            # Open and read the PDB file
            pdb_file_path = os.path.join(protein_folder, protein_name + '.pdb')

            # If a user uploaded pdb file then we prioritise this structure
            if os.path.exists(pdb_file_path):
                with open(pdb_file_path, 'r') as pdb_file:
                    for line in pdb_file:
                        pdb_content += line
            else:
                pred_pdb_file_path = os.path.join(result_folder, protein_name + "_pred_protein.pdb")
                with open(pred_pdb_file_path, 'r') as pdb_file:
                    for line in pdb_file:
                        pdb_content += line

            ligand_names = self.get_ligand_names(result_folder)

            for ligand_name in ligand_names:

                ligand_file = f'{result_folder}/{ligand_name}_pred_ligand.sdf'

                plddt_df = pd.read_csv(f"{result_folder}/{ligand_name}_ligand_plddt.csv", header=None)
                ligand_plddt = np.round(plddt_df[0].mean(),1)

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
                    'affinity': ligand_plddt
                    }
                )

        return results

    def check_result_available(self, experiments_folder, experiment_id: str, protein_id: str, ligand_id: str):
        experiment_folder = os.path.join(experiments_folder, experiment_id)

        protein_folder = os.path.join(experiment_folder, 'results', protein_id)
        result_folder = os.path.join(protein_folder, ligand_id, 'result')

        ligand_file = f'{result_folder}/{ligand_id}_pred_ligand.sdf'

        plddt_df_file = f"{result_folder}/{ligand_id}_ligand_plddt.csv"

        if not (os.path.exists(ligand_file) and os.path.exists(plddt_df_file)):
            return False
        return True

    def get_protein_ids(self, experiments_folder, experiment_id):
        experiment_folder = os.path.join(experiments_folder, experiment_id, 'results')
        if not os.path.exists(experiment_folder):
            os.mkdir(experiment_folder)
        protein_ids = [d for d in os.listdir(experiment_folder) \
        if os.path.isdir(os.path.join(experiment_folder, d))]

        return protein_ids


    def get_ligand_names(self, result_folder):
        ligand_ids = [d for d in os.listdir(result_folder) 
                      if os.path.isdir(os.path.join(result_folder, d))]
        return ligand_ids
