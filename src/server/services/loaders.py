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

class DTILoader:
    def __init__(self):
        pass

    def get_dti_results(self, experiments_folder: str, experiment_id: str):
        experiment_folder = os.path.join(experiments_folder, experiment_id)
        import glob
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

    def get_ligand_names(self, result_folder):
        suffix = '_pred_ligand'
        ligands_path = glob.glob(result_folder + f'/*{suffix}.sdf')
        file_names = [
            Path(ligand_path).stem for ligand_path in ligands_path
        ]
        return [file_name.replace(suffix, '') for file_name in file_names]
