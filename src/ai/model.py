import shutil

import torch
from torch import Tensor
from transformers import AutoModelForSequenceClassification, AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37
from torch_geometric.loader import DataLoader

import numpy as np
import requests
from tqdm import tqdm

from src.ai.exceptions.model_not_loaded_ex import ModelNotLoadedException
from src.ai.custom_models.custom_models import SimpleGOMultiLayerPerceptron, SimpleSolubilityMultiLayerPerceptron
from src.ai.custom_models.drug_target.utils import install_p2rank, read_sdf_files

from Bio.PDB import PDBParser
import rdkit.Chem as Chem
import pandas as pd

from src.ai.custom_models.drug_target.tankbind.feature_utils import get_protein_feature, \
    generate_sdf_from_smiles_using_rdkit, extract_torchdrug_feature_from_mol
from src.ai.custom_models.drug_target.tankbind.generation_utils import get_LAS_distance_constraint_mask, \
    get_info_pred_distance, write_with_new_coords
from src.ai.custom_models.drug_target.tankbind.model import get_model
from src.ai.custom_models.drug_target.tankbind.data import TankBind_prediction

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


class DrugTargetInteraction(BaseModel):
    def __init__(self, model_name, gpu, model_task=""):
        super().__init__(model_name, gpu, model_task)
        if gpu:
            self.device = torch.device('cuda')
        else:
            self.device = torch.device('cpu')

        logger.info("Loading p2rank...")
        install_p2rank()
        logger.info("P2rank has been loaded!")

    def prepare_folders(self, experiment_folder: str, protein_names: List[str]):
        for protein_name in protein_names:
            os.system(f'mkdir -p {experiment_folder}/{protein_name}')

    def prepare_data(
            self,
            ligands_names: List[str],
            ligands_smiles: List[FileStorage],
            protein_names: List[str],
            protein_files: List[FileStorage],
            experiment_folder: str
    ):
        self.prepare_folders(experiment_folder, protein_names)
        self.protein_dict = self._prepare_protein_data(experiment_folder, protein_files, protein_names)
        self.protein2ligandcompdict = {}
        for protein_name in protein_names:
            ligand2compounddict = {}
            for ligand_name, ligand_smiles in zip(ligands_names, ligands_smiles):
                # getting compund feature
                compound_dict = {}
                rdkitMolFile = f"{experiment_folder}/{protein_name}/{ligand_name}_mol_from_rdkit.sdf"
                shift_dis = 0  # for visual only, could be any number, shift the ligand away from the protein.
                generate_sdf_from_smiles_using_rdkit(ligand_smiles, rdkitMolFile, shift_dis=shift_dis)
                mol = Chem.MolFromMolFile(rdkitMolFile)
                compound_dict[protein_name + "_" + f"{ligand_name}_rdkit"] = extract_torchdrug_feature_from_mol(mol,
                                                                                                                has_LAS_mask=True)
                ligand2compounddict[ligand_name] = [compound_dict, rdkitMolFile]
                info = []
                for compound_name in list(compound_dict.keys()):
                    # use protein center as the block center.
                    com = ",".join([str(a.round(3)) for a in self.protein_dict[protein_name][0].mean(axis=0).numpy()])
                    info.append([protein_name, compound_name, "protein_center", com])
                    p2rankFile = f"{experiment_folder}/{protein_name}/p2rank/{protein_name}.pdb_predictions.csv"
                    pocket = pd.read_csv(p2rankFile)
                    pocket.columns = pocket.columns.str.strip()
                    pocket_coms = pocket[['center_x', 'center_y', 'center_z']].values
                    for ith_pocket, com in enumerate(pocket_coms):
                        com = ",".join([str(a.round(3)) for a in com])
                        info.append([protein_name, compound_name, f"pocket_{ith_pocket + 1}", com])
                info = pd.DataFrame(info, columns=['protein_name', 'compound_name', 'pocket_name', 'pocket_com'])
                info.to_csv(f"{experiment_folder}/{protein_name}/{ligand_name}_temp_info.csv")

            self.protein2ligandcompdict[protein_name] = ligand2compounddict

    def _prepare_protein_data(self, experiment_folder, protein_files: List[FileStorage], protein_names):
        protein_dict = {}
        for protein_file, protein_name in zip(protein_files, protein_names):
            # p2rank
            original_protein_file = os.path.join(experiment_folder, protein_file.filename)
            destination_protein_file = os.path.join(experiment_folder, protein_name, protein_file.filename)
            shutil.copyfile(original_protein_file, destination_protein_file)
            ds = f"{experiment_folder}/{protein_name}/protein_list.ds"
            with open(ds, "w") as out:
                out.write(f"{protein_file.filename}\n")
            p2rank_exec = dirname(os.path.abspath(__file__)) + "/custom_models/drug_target/tankbind/p2rank_2.4.1/prank"
            print("P2RANK_EXEC: " + p2rank_exec)
            p2rank = f"bash {p2rank_exec}"
            cmd = f"{p2rank} predict {ds} -o {experiment_folder}/{protein_name}/p2rank -threads 1"
            os.system(cmd)
            # getting protein feature
            parser = PDBParser(QUIET=True)
            s = parser.get_structure("x", os.path.join(experiment_folder, protein_file.filename))
            res_list = list(s.get_residues())
            protein_dict[protein_name] = get_protein_feature(res_list)
        return protein_dict

    def load_model(self):
        logging.basicConfig(level=logging.INFO)
        self.model = get_model(0, logging, self.device)
        modelFile = dirname(os.path.abspath(__file__)) + "/custom_models/drug_target/models_ckpt/self_dock.pt"
        self.model.load_state_dict(torch.load(modelFile, map_location=self.device))
        self.model.eval()

    def _raw_inference(self, experiment_folder, protein_names, ligands_names):

        def post_process(chosen, rdkitMolFile, dataset, result_folder):
            for i, line in chosen.iterrows():
                idx = line['index']
                pocket_name = line['pocket_name']
                compound_name = line['compound_name']
                ligandName = compound_name.split("_")[1]
                coords = dataset[idx].coords.to(self.device)
                protein_nodes_xyz = dataset[idx].node_xyz.to(self.device)
                n_compound = coords.shape[0]
                n_protein = protein_nodes_xyz.shape[0]
                y_pred = self.y_pred_list[idx].reshape(n_protein, n_compound).to(self.device)
                y = dataset[idx].dis_map.reshape(n_protein, n_compound).to(self.device)
                compound_pair_dis_constraint = torch.cdist(coords, coords)
                mol = Chem.MolFromMolFile(rdkitMolFile)
                LAS_distance_constraint_mask = get_LAS_distance_constraint_mask(mol).bool()
                info = get_info_pred_distance(coords, y_pred, protein_nodes_xyz, compound_pair_dis_constraint,
                                              LAS_distance_constraint_mask=LAS_distance_constraint_mask,
                                              n_repeat=1, show_progress=False)
                toFile = f'{result_folder}/{ligandName}_tankbind.sdf'
                new_coords = info.sort_values("loss")['coords'].iloc[0].astype(np.double)
                write_with_new_coords(mol, new_coords, toFile)

        for protein_idx in range(len(protein_names)):
            protein_name = protein_names[protein_idx]
            logger.info(f"Making dti predictions for {protein_idx+1} out of {len(protein_names)} protein...")
            result_folder = os.path.join(experiment_folder, protein_name, 'result')
            if not os.path.exists(result_folder):
                os.mkdir(result_folder)
            ligand2compounddict = self.protein2ligandcompdict[protein_name]
            for ligand_idx in range(len(ligands_names)):
                ligand_name = ligands_names[ligand_idx]
                logger.info(f"Making dti predictions for {protein_idx+1} out of {len(protein_names)} protein...\n \
                Currently predicting {ligand_idx+1} out of {len(ligands_names)} ligands...")
                info = pd.read_csv(os.path.join(experiment_folder,protein_name, f"{ligand_name}_temp_info.csv"))
                compound_dict, rdkitMolFile = ligand2compounddict[ligand_name]
                dataset_path = f"{experiment_folder}/{protein_name}/{ligand_name}_dataset/"
                os.system(f"rm -r {dataset_path}")
                os.system(f"mkdir -p {dataset_path}")
                dataset = TankBind_prediction(dataset_path, data=info, protein_dict=self.protein_dict,
                                              compound_dict=compound_dict)
                batch_size = 1
                data_loader = DataLoader(dataset, batch_size=batch_size, follow_batch=['x', 'y', 'compound_pair'],
                                         shuffle=False, num_workers=8)
                affinity_pred_list = []
                self.y_pred_list = []
                for data in tqdm(data_loader):
                    data = data.to(self.device)
                    y_pred, affinity_pred = self.model(data)
                    affinity_pred_list.append(affinity_pred.detach().cpu())
                    for i in range(data.y_batch.max() + 1):
                        self.y_pred_list.append((y_pred[data['y_batch'] == i]).detach().cpu())
                affinity_pred_list = torch.cat(affinity_pred_list)
                info['affinity'] = affinity_pred_list
                info.to_csv(f"{result_folder}/{ligand_name}_info_with_predicted_affinity.csv")
                chosen = info.loc[
                    info.groupby(['protein_name', 'compound_name'], sort=False)['affinity'].agg('idxmax')].reset_index()
                post_process(chosen, rdkitMolFile, dataset=dataset, result_folder=result_folder)

    def predict(self, ligand_files: List[FileStorage], protein_files: List[FileStorage], experiment_folder: str):
        logger.info("Making dti predictions...")
        ligand_files = [os.path.join(experiment_folder, file.filename) for file in ligand_files]
        ligands_names, ligands_smiles = read_sdf_files(ligand_files)
        protein_names = [os.path.splitext(x.filename)[0] for x in protein_files]

        self.prepare_data(ligands_names, ligands_smiles, protein_names, protein_files, experiment_folder)
        self._raw_inference(experiment_folder, protein_names, ligands_names)
        logger.info("DTI predictions are completed!")
