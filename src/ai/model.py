import shutil
import subprocess
import uuid

import torch
from torch import Tensor
from transformers import AutoModelForSequenceClassification, AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37
from torch_geometric.loader import DataLoader

import numpy as np
import requests
from tqdm import tqdm
import pickle
from pathlib import Path

from src.server import settings
from src.server.services.progress import ProgressTracker
from src.ai.exceptions.model_not_loaded_ex import ModelNotLoadedException
from src.ai.custom_models.custom_models import SimpleGOMultiLayerPerceptron, SimpleSolubilityMultiLayerPerceptron
from src.ai.custom_models.drug_target.utils import install_p2rank, read_sdf_files, get_sequence

from src.ai.custom_models.drug_target.umol.src.check_msa_colab import process_a3m
from src.ai.custom_models.drug_target.umol.src.make_ligand_feats_colab import bonds_from_smiles
from src.ai.custom_models.drug_target.umol.src.make_msa_seq_feats_colab import process
from src.ai.custom_models.drug_target.umol.src.net.model import config
from src.ai.custom_models.drug_target.umol.src.predict_colab import predict
from src.ai.custom_models.drug_target.umol.src.relax.align_ligand_conformer_colab import read_pdb, \
    generate_best_conformer, align_coords_transform, write_sdf

from Bio.PDB import PDBParser
import rdkit.Chem as Chem
import pandas as pd

from esm import FastaBatchedDataset, pretrained

from typing import List, Tuple
import os
import logging
from src.server.initializers.loggers import logger
from argparse import Namespace
from werkzeug.datastructures import FileStorage

dirname = os.path.dirname


class BaseModel:
    def __init__(self, model_name: str, gpu: bool, model_task=""):
        self.model_name = model_name
        self.model_task = model_task
        self.model = None
        self.tokenizer = None
        self.gpu = gpu

    def load_model(self):
        """Load model and tokenizer here"""
        pass

    # Method to get raw model outputs
    def _raw_inference(self, input: str):
        pass

    # Method to return raw outputs in the desired format
    def predict(self, input: str):
        pass


class ClassificationModel(BaseModel):
    def __init__(self, model_name: str, gpu: bool, model_task=""):
        super().__init__(model_name, gpu, model_task)

    def set_labels(self, labels: List[str]):
        self.labels = labels

    def load_model(self):
        logger.info("Loading classification model")
        self.model = AutoModelForSequenceClassification.from_pretrained(self.model_name)
        self.tokenizer = AutoTokenizer.from_pretrained(self.model_name)
        self.model.eval()
        return super().load_model()

    def _raw_inference(self, sequence: str) -> Tensor:
        inputs = self.tokenizer(sequence, return_tensors='pt', padding=True, truncation=True)
        with torch.no_grad():
            outputs = self.model(**inputs)
            probabilities = torch.nn.functional.softmax(outputs.logits, dim=-1).tolist()[0]
        return probabilities

    def predict(self, sequence: str) -> List[List[Tuple[str, int]]]:
        logger.info("Making classification predictions...")
        if not self.tokenizer or not self.model:
            raise ModelNotLoadedException()

        probabilities = self._raw_inference(sequence)
        prob_table = {key: value for key, value in list(zip(self.labels, probabilities))}
        logger.info("Successfully made classification predictions!")
        return prob_table


class Folding(BaseModel):

    def __init__(self, model_name: str, gpu: bool, model_task=""):
        super().__init__(model_name, gpu, model_task)
        self.model_name = model_name
        self.gpu = gpu

    def load_model(self):
        logger.info("Loading folding model")
        self.tokenizer = AutoTokenizer.from_pretrained(self.model_name)
        self.model = EsmForProteinFolding.from_pretrained(self.model_name, low_cpu_mem_usage=True)
        if self.gpu:
            self.model = self.model.cuda()
        self.model.eval()

        return super().load_model()

    def convert_outputs_to_pdb(self, outputs) -> List[str]:
        final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)
        outputs = {k: v.to("cpu").numpy() for k, v in outputs.items()}
        final_atom_positions = final_atom_positions.cpu().numpy()
        final_atom_mask = outputs["atom37_atom_exists"]
        pdbs = []
        for i in range(outputs["aatype"].shape[0]):
            aa = outputs["aatype"][i]
            pred_pos = final_atom_positions[i]
            mask = final_atom_mask[i]
            resid = outputs["residue_index"][i] + 1
            pred = OFProtein(
                aatype=aa,
                atom_positions=pred_pos,
                atom_mask=mask,
                residue_index=resid,
                b_factors=outputs["plddt"][i],
                chain_index=outputs["chain_index"][i] if "chain_index" in outputs else None,
            )
            pdbs.append(to_pdb(pred))
        return pdbs

    def _raw_inference(self, sequence: str):

        tokenized_input = self.tokenizer([sequence], return_tensors="pt", add_special_tokens=False)['input_ids']

        if self.gpu:
            tokenized_input = tokenized_input.cuda()

        output = None

        with torch.no_grad():
            output = self.model(tokenized_input)

        return output

    def predict(self, sequence: str) -> str:
        logger.info("Predicting folded structure...")
        if not self.tokenizer or not self.model:
            raise ModelNotLoadedException()
        output = self._raw_inference(sequence)

        pdbs = self.convert_outputs_to_pdb(output)

        logger.info("Successfully predicted folded structure!")
        return "".join(pdbs)


class ESM2EmbeddingGenerator(BaseModel):
    def __init__(self, model_name, gpu, model_task=""):
        super().__init__(model_name, gpu)
        if gpu:
            self.device = torch.device('cuda')
        else:
            self.device = torch.device('cpu')

    def load_model(self):
        self.model, self.alphabet = pretrained.load_model_and_alphabet(self.model_name)
        self.model.to(self.device)
        self.model.eval()

    def predict(self, protein_sequence: str):
        include = "mean"
        repr_layers = [-1]
        truncation_seq_length = 1022
        nogpu = False
        toks_per_batch = 4096

        args = Namespace(
            repr_layers=repr_layers,
            include=include,
            truncation_seq_length=truncation_seq_length,
            nogpu=nogpu,
            toks_per_batch=toks_per_batch
        )

        if self.device == torch.device('cuda'):
            self.model = self.model.cuda()
            print("Transferred model to GPU")

        dataset = FastaBatchedDataset(["mockId"], [protein_sequence])
        batches = dataset.get_batch_indices(args.toks_per_batch, extra_toks_per_seq=1)
        data_loader = torch.utils.data.DataLoader(
            dataset, collate_fn=self.alphabet.get_batch_converter(args.truncation_seq_length), batch_sampler=batches
        )

        return_contacts = "contacts" in args.include

        repr_layers = [(i + self.model.num_layers + 1) % (self.model.num_layers + 1) for i in args.repr_layers]

        with torch.no_grad():
            for batch_idx, (labels, strs, toks) in enumerate(data_loader):

                if self.device == torch.device('cuda'):
                    toks = toks.to(device="cuda", non_blocking=True)

                out = self.model(toks, repr_layers=repr_layers, return_contacts=return_contacts)

                logits = out["logits"].to(device="cpu")
                representations = {
                    layer: t.to(device="cpu") for layer, t in out["representations"].items()
                }

                for i, label in enumerate(labels):
                    label = label.split()[0]
                    result = {"label": label}
                    truncate_len = min(args.truncation_seq_length, len(strs[i]))
                    result["mean_representations"] = {
                        layer: t[i, 1: truncate_len + 1].mean(0).clone() \
                        for layer, t in representations.items()
                    }
                    return result["mean_representations"][30]


class GeneOntologyPrediction(BaseModel):
    """
    Look for a better model there https://www.kaggle.com/competitions/cafa-5-protein-function-prediction
    """

    def __init__(self, model_name, gpu, model_task=""):
        super().__init__(model_name, gpu, model_task)
        if gpu:
            self.device = torch.device('cuda')
        else:
            self.device = torch.device('cpu')
        self.embedding_model = None

    def load_model(self):
        logger.info("Loading gene ontology model...")
        self.model = SimpleGOMultiLayerPerceptron(input_dim=640, num_classes=200).to(self.device)
        self._load_model_state_dict()
        self.model.eval()
        self.labels = np.load(dirname(os.path.abspath(__file__)) + \
                              "/custom_models/models/gene_ontology/go_labels_200.npy", allow_pickle=True)
        logger.info("Gene ontology model has been loaded!")

    def set_embedding_model(self, embedding_model):
        self.embedding_model = embedding_model

    def _load_model_state_dict(self):
        # URL of the model's .pth file on Hugging Face Model Hub
        model_url = "https://huggingface.co/thomasshelby/go_prediction/resolve/main/go_model_150M.pth"

        # Path where you want to save the downloaded .pth file
        local_path = dirname(os.path.abspath(__file__)) + "/custom_models/models/gene_ontology/go_model_150M.pth"

        if not os.path.exists(local_path):
            with open(local_path, 'wb') as f:
                response = requests.get(model_url)
                f.write(response.content)

        self.model.load_state_dict(torch.load(local_path, map_location=self.device))

    def predict(self, protein_id: str):
        logger.info("Making gene ontology predictions...")
        embedding = self.embedding_model.predict(protein_id)
        embedding = embedding.to(self.device)
        raw_outputs = torch.nn.functional.sigmoid(self.model(embedding)).squeeze().detach().cpu().numpy()

        outputs = {label: confidence for label, confidence in zip(self.labels, raw_outputs)}
        logger.info("Successfully predicted gene ontology!")
        return outputs


class SolubilityPrediction(BaseModel):
    def __init__(self, model_name, gpu, model_task=""):
        super().__init__(model_name, gpu, model_task)
        if gpu:
            self.device = torch.device('cuda')
        else:
            self.device = torch.device('cpu')
        self.embedding_model = None

    def load_model(self):
        logger.info("Loading solubility model...")
        self.model = SimpleSolubilityMultiLayerPerceptron().to(self.device)
        self._load_model_state_dict()
        self.model.eval()
        logger.info("Successfully loaded solubility model!")

    def set_embedding_model(self, embedding_model):
        self.embedding_model = embedding_model

    def _load_model_state_dict(self):
        # URL of the model's .pth file on Hugging Face Model Hub
        model_url = "https://huggingface.co/thomasshelby/solubility_model/resolve/main/solubility_model.pth"

        # Path where you want to save the downloaded .pth file
        local_path = dirname(os.path.abspath(__file__)) + "/custom_models/models/solubility/solubility_model.pth"

        if not os.path.exists(local_path):
            with open(local_path, 'wb') as f:
                response = requests.get(model_url)
                f.write(response.content)

        self.model.load_state_dict(torch.load(local_path, map_location=self.device))

    def predict(self, protein_id: str):
        logger.info("Predicting solubility...")
        embedding = self.embedding_model.predict(protein_id)
        embedding = embedding.to(self.device)
        outputs = self.model(embedding.float())
        logger.info("Solubility has been predicted!")
        return {'solubility': outputs.item()}




class PocketPredictor(BaseModel):
    def __init__(self, model_name: str, model_task=""):
        self.model_name = model_name
        self.model_task = model_task

    def load_model(self):
        """Load model and tokenizer here"""
        install_p2rank()

    # Method to get raw model outputs
    def _raw_inference(self, protein_file_path: str, save_dir: str):
        protein_filename = os.path.split(protein_file_path)[1]
        ds = f"{save_dir}/protein_list.ds"
        with open(ds, "w") as out:
            out.write(f"{protein_filename}\n")
        p2rank_exec = os.path.join(os.path.dirname(os.path.abspath(__file__)), "custom_models/drug_target/p2rank_2.4.1/prank")
        cmd = ["bash", p2rank_exec, "predict", ds, "-o", f"{save_dir}/p2rank", "-threads", "1"]

        # Run the command and wait for it to complete
        subprocess.run(cmd, check=True)

        p2rankFile = os.path.join(save_dir, 'p2rank', f"{Path(protein_filename).stem}.pdb_predictions.csv")
        pocket = pd.read_csv(p2rankFile, skipinitialspace=True)

        if pocket.empty:
            return np.array([])

        residue_ids = pocket['residue_ids'].str.split()
        all_ids = [int(item.split('_')[1]) for sublist in residue_ids for item in sublist]

        pocket_ids = np.sort(np.array(all_ids))

        np.save(os.path.join(save_dir, "pocket.npy"), pocket_ids)

        pocket_ids = pocket_ids.tolist()

        return pocket_ids

        
    # Method to return raw outputs in the desired format
    def predict(self, protein_file_path: str, save_dir: str):
        return self._raw_inference(protein_file_path, save_dir)



class DrugTargetInteraction(BaseModel):
    def __init__(self, model_name, gpu, model_task=""):
        super().__init__(model_name, gpu, model_task)

    def submit_fasta_and_save_a3m(self, api_url, fasta_file_path, save_dir):
        # Extract the base name to create the .a3m file name
        base_name = os.path.splitext(os.path.basename(fasta_file_path))[-2]
        a3m_file_path = os.path.join(save_dir, f"{base_name}.a3m")

        # Open the FASTA file and send it to the API
        with open(fasta_file_path, 'rb') as file:
            files = {'sequence_file': file}
            response = requests.post(api_url, files=files)

        # Check if the request was successful
        if response.status_code == 200:
            # Write the .a3m content to a file
            with open(a3m_file_path, 'w') as a3m_file:
                a3m_file.write(response.json()['alignment'])
            logger.info(f"Saved .a3m file for {fasta_file_path} as {a3m_file_path}")
        else:
            logger.error(f"Failed to get .a3m for {fasta_file_path}: {response.status_code}")


    def prepare_ligand_data(self, fasta_file, ligands, ligand_names, save_dir):
        for ligand, ligand_name in zip(ligands, ligand_names):
            protein_name = os.path.splitext(os.path.basename(fasta_file))[0]

            if not os.path.exists(save_dir):
                os.makedirs(save_dir)

            MSA = os.path.join(save_dir, protein_name + '.a3m')
            PROCESSED_MSA = os.path.join(save_dir, protein_name + '_processed.a3m')
            process_a3m(MSA, get_sequence(fasta_file), PROCESSED_MSA)
            MSA = PROCESSED_MSA

            # Process MSA features
            feature_dict = process(fasta_file, [MSA])  # Assuming MSA is defined elsewhere
            features_output_path = os.path.join(save_dir, 'msa_features.pkl')
            with open(features_output_path, 'wb') as f:
                pickle.dump(feature_dict, f, protocol=4)
            logger.info('Saved MSA features to', features_output_path)

            atom_encoding = {'B': 0, 'C': 1, 'F': 2, 'I': 3, 'N': 4, 'O': 5, 'P': 6, 'S': 7, 'Br': 8, 'Cl': 9,
                             # Individual encoding
                             'As': 10, 'Co': 10, 'Fe': 10, 'Mg': 10, 'Pt': 10, 'Rh': 10, 'Ru': 10, 'Se': 10, 'Si': 10,
                             'Te': 10, 'V': 10, 'Zn': 10
                             # Joint (rare)
                             }

            # Process ligand features
            atom_types, atoms, bond_types, bond_lengths, bond_mask = bonds_from_smiles(ligand, atom_encoding)
            ligand_inp_feats = {
                'atoms': atoms,
                'atom_types': atom_types,
                'bond_types': bond_types,
                'bond_lengths': bond_lengths,
                'bond_mask': bond_mask
            }

            features_output_path = os.path.join(save_dir, f'{ligand_name}_ligand_inp_features.pkl')
            with open(features_output_path, 'wb') as f:
                pickle.dump(ligand_inp_feats, f, protocol=4)
            print('Saved ligand features to', features_output_path)

    def load_model(self):
        print("Loading uMol DTI params...")
        self.model_params_path = os.path.join(dirname(os.path.abspath(__file__)), 'custom_models', 'drug_target',
                                              'models_ckpt', 'params.pkl')
        url = 'https://huggingface.co/thomasshelby/uMol_params/resolve/main/params.pkl?download=true'
        if not os.path.exists(self.model_params_path):
            response = requests.get(url, stream=True)
            if response.status_code == 200:
                with open(self.model_params_path, 'wb') as f:
                    f.write(response.content)
            else:
                raise Exception(f"Failed to download file from {url}")

        pass

    def _raw_inference(self, protein_file_paths, protein_names, protein_ids, ligands, ligands_names, experiment_folder, binding_pockets, num_recycles):
        results_dir = os.path.join(experiment_folder, 'results')
        if not os.path.exists(results_dir):
                os.makedirs(results_dir)
        experiment_progress_tracker = ProgressTracker(experiment_folder, protein_names)
        for protein_file_path, protein_name, protein_id, binding_pocket in zip(protein_file_paths, protein_names, protein_ids, binding_pockets):
            result_id = protein_id
            save_dir = os.path.join(results_dir, result_id)
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            #self.submit_fasta_and_save_a3m(settings.FASTA_API, protein_file_path, save_dir)
            self.prepare_ligand_data(protein_file_path, ligands, ligands_names, save_dir)

            for LIGAND, LIGAND_NAME in zip(ligands, ligands_names):
                protein_folder = save_dir
                result_folder = os.path.join(protein_folder, 'result')

                protein_progress_tracker = ProgressTracker(protein_folder, tasks=ligands_names)

                if not os.path.exists(result_folder):
                    os.makedirs(result_folder)

                result_folder = os.path.join(protein_folder, 'result/')

                MSA_FEATS = os.path.join(protein_folder, 'msa_features.pkl')
                LIGAND_FEATS = os.path.join(protein_folder, f'{LIGAND_NAME}_ligand_inp_features.pkl')

                with open(self.model_params_path, 'rb') as file:
                    PARAMS = pickle.load(file)

                protein_name = protein_name.split('/')[-1]
                # Predict
                predict(config.CONFIG, MSA_FEATS, LIGAND_FEATS, protein_name, binding_pocket, PARAMS, num_recycles, outdir=result_folder)

                # Process the prediction
                RAW_PDB = os.path.join(result_folder, f'{protein_name}_pred_raw.pdb')

                # Get a conformer
                pred_ligand = read_pdb(RAW_PDB)
                best_conf, best_conf_pos, best_conf_err, atoms, nonH_inds, mol, best_conf_id = generate_best_conformer(pred_ligand['chain_coords'], LIGAND)

                # Align it to the prediction
                aligned_conf_pos = align_coords_transform(pred_ligand['chain_coords'], best_conf_pos, nonH_inds)

                # Write sdf
                sdf_output_path = os.path.join(result_folder, f'{LIGAND_NAME}_pred_ligand.sdf')
                write_sdf(mol, best_conf, aligned_conf_pos, best_conf_id, sdf_output_path)

                # Extract ATOM and HETATM records from the PDB file
                protein_pdb_path = os.path.join(result_folder, f'{protein_name}_pred_protein.pdb')
                ligand_plddt_path = os.path.join(result_folder, f'{LIGAND_NAME}_ligand_plddt.csv')

                with open(RAW_PDB, 'r') as infile, open(protein_pdb_path, 'w') as protein_out, open(ligand_plddt_path, 'w') as ligand_out:
                    for line in infile:
                        if line.startswith('ATOM'):
                            protein_out.write(line)
                        elif line.startswith('HETATM'):
                            ligand_out.write(line[64:66] + '\n')  # Extracting plDDT values

                #protein_progress_tracker.update_progress(LIGAND_NAME)

            #experiment_progress_tracker.update_progress(protein_name)


    def predict(self, ligand_files_paths, protein_files, protein_ids, binding_pockets, experiment_folder: str):
        logger.info("Making dti predictions...")
        ligands_names, ligands_smiles = read_sdf_files(ligand_files_paths)
        protein_file_paths = [os.path.splitext(x)[0] for x in protein_files]

        protein_names = [os.path.basename(protein_file_path) for protein_file_path in protein_file_paths]

        self._raw_inference(
            protein_file_paths=protein_files,
            protein_names=protein_names,
            protein_ids=protein_ids,
            ligands=ligands_smiles,
            ligands_names=ligands_names,
            experiment_folder=experiment_folder,
            binding_pockets=binding_pockets,
            num_recycles=3)
