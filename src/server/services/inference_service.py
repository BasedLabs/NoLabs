import json
import os
from typing import Dict, List

from rdkit import Chem
from rdkit.Chem import SDMolSupplier

from src.ai.model_factory import create_model
from src.ai.pipeline import Pipeline

def create_pipeline(use_gpu=False, is_test = False) -> Pipeline:
    pipeline = Pipeline()
    models_metadata = get_models_from_config(is_test)

    for model_metadata in models_metadata:
        model = create_model(model_metadata, use_gpu)
        pipeline.add_model(model)

    return pipeline

def get_localisation_output(pipeline, amino_acid_sequence: str) -> Dict:
    assert amino_acid_sequence
    model = pipeline.get_model_by_task("localisation")
    return model.predict(amino_acid_sequence)

def get_folding_output(pipeline, amino_acid_sequence: str) -> str:
    assert amino_acid_sequence
    model = pipeline.get_model_by_task("folding")
    return model.predict(amino_acid_sequence)

def get_gene_ontology_output(pipeline, amino_acid_sequence: str) -> Dict:
    model = pipeline.get_model_by_task("gene_ontology")
    res = model.predict(amino_acid_sequence)
    return res

def get_solubility_output(pipeline, amino_acid_sequence: str) -> Dict:
    model = pipeline.get_model_by_task("solubility")
    res = model.predict(amino_acid_sequence)
    return res

def generate_dti_results(pipeline, ligand_files: List[str], protein_file: str):
    model = pipeline.get_model_by_task("dti")
    model.predict(ligand_files, protein_file)

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

    ligand_names = [os.path.splitext(file.filename)[0] for file in ligand_files]

    for ligand_name in ligand_names:

        ligand_file = f'{result_folder}/{ligand_name}_tankbind.sdf'


        sdf_supplier = SDMolSupplier(ligand_file)
        # Initialize an empty string to store the SDF contents
        sdf_contents = ""
        # Iterate through the molecules in the SDF file and append their representations to the string
        for mol in sdf_supplier:
            if mol is not None:
                # Convert the molecule to an SDF block and append it to the string
                sdf_contents += Chem.MolToMolBlock(mol) + "\n"
        ligands_sdf_contents.append(sdf_contents)

    return pdb_content, protein_name, ligands_sdf_contents, ligand_names

def save_uploaded_files(pipeline, files):
    model = pipeline.get_model_by_task("dti")
    save_folder = model.experiment_folder

    # Create the 'temp/' folder if it doesn't exist
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)

    saved_file_paths = []

    for file in files:
        # Generate a unique filename to avoid overwriting
        filename = os.path.join(save_folder, file.filename)

        # Save the file to the 'temp/' folder
        file.save(filename)

        saved_file_paths.append(filename)

    return saved_file_paths

def get_pipeline_output(pipeline, amino_acid_sequence: str) -> str:
    assert amino_acid_sequence
    return pipeline.predict(amino_acid_sequence)

def read_json_file(file_path):
    with open(file_path, 'r') as json_file:
        data = json.load(json_file)
    return data

def read_config():
    dirname = os.path.dirname
    root_directory = dirname(dirname(dirname(dirname(os.path.abspath(__file__)))))
    file_path = os.path.join(root_directory, 'config.json')
    data = read_json_file(file_path)
    return data

def read_test_config():
    dirname = os.path.dirname
    root_directory = dirname(dirname(dirname(dirname(os.path.abspath(__file__)))))
    file_path = os.path.join(root_directory, 'test/ai/mock_resources/config.json')
    data = read_json_file(file_path)
    return data

def get_models_from_config(is_test):
    data = None
    if is_test:
        data = read_test_config()
    else:
        data = read_config()
    models = data.get('models', [])
    return models
