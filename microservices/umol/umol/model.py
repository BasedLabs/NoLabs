import shutil
import sys
import os
from os.path import dirname

sources_root = os.path.join(dirname(__file__), 'umol_source')
print(f'Additional sources root {sources_root}')
sys.path.append(sources_root)

from typing import List, Union, Tuple
from fastapi import UploadFile

import pickle

from umol.umol_source.src.check_msa_colab import process_a3m
from umol.umol_source.src.make_ligand_feats_colab import bonds_from_smiles
from umol.umol_source.src.make_msa_seq_feats_colab import process
from umol.umol_source.src.net.model import config
from umol.umol_source.src.predict_colab import predict
from umol.umol_source.src.relax.align_ligand_conformer_colab import read_pdb, \
    generate_best_conformer, align_coords_transform, write_sdf


class DrugTargetInteraction:
    def __init__(self):
        self.model_params_path = ''
        self.load_model()

    def prepare_ligand_data(self, ligand, save_dir):
        atom_encoding = {'B': 0, 'C': 1, 'F': 2, 'I': 3, 'N': 4, 'O': 5, 'P': 6, 'S': 7, 'Br': 8, 'Cl': 9,
                         # Individual encoding
                         'As': 10, 'Co': 10, 'Fe': 10, 'Mg': 10, 'Pt': 10, 'Rh': 10, 'Ru': 10, 'Se': 10, 'Si': 10,
                         'Te': 10, 'V': 10, 'Zn': 10
                         # Joint (rare)
                         }

        # Process ligand modules
        atom_types, atoms, bond_types, bond_lengths, bond_mask = bonds_from_smiles(ligand, atom_encoding)
        ligand_inp_feats = {
            'atoms': atoms,
            'atom_types': atom_types,
            'bond_types': bond_types,
            'bond_lengths': bond_lengths,
            'bond_mask': bond_mask
        }

        features_output_path = os.path.join(save_dir, 'ligand_inp_features.pkl')
        with open(features_output_path, 'wb') as f:
            pickle.dump(ligand_inp_feats, f, protocol=4)
        print('Saved ligand modules to', features_output_path)

    def load_model(self):
        print("Loading uMol DTI params...")
        self.model_params_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'umol_source', 'params.pkl')

    def _raw_inference(self,
                       protein_fasta_path: str,
                       ligand_smiles: str,
                       msa_file_path: str,
                       binding_pocket: List[int],
                       save_dir: Union[str, bytes]) -> Tuple[str, str, List[int]]:

        MSA = msa_file_path
        PROCESSED_MSA = os.path.join(save_dir, 'protein_processed.a3m')
        process_a3m(MSA, get_sequence(protein_fasta_path), PROCESSED_MSA)
        MSA = PROCESSED_MSA

        # Process MSA modules
        feature_dict = process(protein_fasta_path, [MSA])  # Assuming MSA is defined elsewhere
        features_output_path = os.path.join(save_dir, 'msa_features.pkl')
        with open(features_output_path, 'wb') as f:
            pickle.dump(feature_dict, f, protocol=4)

        self.prepare_ligand_data(ligand_smiles, save_dir)

        result_folder = os.path.join(save_dir, 'result/')

        if not os.path.exists(result_folder):
            os.makedirs(result_folder)

        MSA_FEATS = os.path.join(save_dir, 'msa_features.pkl')
        LIGAND_FEATS = os.path.join(save_dir, 'ligand_inp_features.pkl')

        with open(self.model_params_path, 'rb') as file:
            PARAMS = pickle.load(file)

        # Predict
        predict(config.CONFIG, MSA_FEATS, LIGAND_FEATS, 'protein', binding_pocket, PARAMS, 3,
                outdir=result_folder)

        # Process the prediction
        RAW_PDB = os.path.join(result_folder, f'protein_pred_raw.pdb')

        # Get a conformer
        pred_ligand = read_pdb(RAW_PDB)
        best_conf, best_conf_pos, best_conf_err, atoms, nonH_inds, mol, best_conf_id = generate_best_conformer(
            pred_ligand['chain_coords'], ligand_smiles)

        # Align it to the prediction
        aligned_conf_pos = align_coords_transform(pred_ligand['chain_coords'], best_conf_pos, nonH_inds)

        # Write sdf
        sdf_output_path = os.path.join(result_folder, f'pred_ligand.sdf')
        write_sdf(mol, best_conf, aligned_conf_pos, best_conf_id, sdf_output_path)

        # Extract ATOM and HETATM records from the PDB file
        protein_pdb_path = os.path.join(result_folder, f'pred_protein.pdb')
        ligand_plddt_path = os.path.join(result_folder, f'ligand_plddt.csv')

        with open(RAW_PDB, 'r') as infile, open(protein_pdb_path, 'w') as protein_out, open(ligand_plddt_path,
                                                                                            'w') as ligand_out:
            for line in infile:
                if line.startswith('ATOM'):
                    protein_out.write(line)
                elif line.startswith('HETATM'):
                    ligand_out.write(line[64:66] + '\n')  # Extracting plDDT values

        # Read and return the contents of pred_ligand.sdf
        with open(sdf_output_path, 'r') as f:
            sdf_contents = f.read()

        with open(protein_pdb_path, 'r') as f:
            pdb_contents = f.read()

        # Read and return the pLDDT values from ligand_plddt.csv
        plddt_values = []
        with open(ligand_plddt_path, 'r') as f:
            for line in f:
                try:
                    plddt_values.append(int(line.strip()))
                except ValueError:
                    pass  # Skip lines that can't be converted to integer

        return sdf_contents, pdb_contents, plddt_values

    def predict(self, protein_sequence: str,
                ligand_smiles: str,
                msa_content: str,
                binding_pocket: List[int]) -> Tuple[str, str, List[int]]:

        temp_directory = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'temp')

        if not os.path.exists(temp_directory):
            os.makedirs(temp_directory)

        # Define paths for the FASTA and MSA files
        temp_protein_fasta_path = os.path.join(temp_directory, 'protein.fasta')
        temp_msa_file_path = os.path.join(temp_directory, 'msa.a3m')

        # Write the protein sequence and MSA content to their respective files
        with open(temp_protein_fasta_path, 'w') as fasta_file:
            fasta_file.write(protein_sequence)

        with open(temp_msa_file_path, 'w') as msa_file:
            msa_file.write(msa_content)

        # Get the results from _raw_inference
        sdf_contents, pdb_contents, plddt_values = self._raw_inference(
            protein_fasta_path=temp_protein_fasta_path,
            ligand_smiles=ligand_smiles,
            msa_file_path=temp_msa_file_path,
            binding_pocket=binding_pocket,
            save_dir=temp_directory)

        return sdf_contents, pdb_contents, plddt_values


def get_sequence(file):
    sequence = ''
    with open(file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                continue  # Skip the header line
            sequence += line.strip()  # Append the sequence removing any whitespace
    return sequence
