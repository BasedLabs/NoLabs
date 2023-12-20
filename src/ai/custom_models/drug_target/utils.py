import os
import tarfile
import urllib.request
import rdkit.Chem as Chem
from rdkit.Chem import SDMolSupplier

dirname = os.path.dirname

def read_sdf_files(ligand_files_paths):
    # Initialize lists to store file names and SMILES strings
    file_names = []
    smiles_list = []

    # Iterate through files in the 'temp/' folder
    for file_path in ligand_files_paths:

        # Extract file name without the '.sdf' extension
        file_name_without_extension = os.path.splitext(os.path.basename(file_path))[0]

        # Initialize an SDMolSupplier to read the SDF file
        sdf_supplier = SDMolSupplier(file_path)

        # Iterate through molecules in the SDF file and convert to SMILES
        for mol in sdf_supplier:
            if mol is not None:
                # Convert the molecule to SMILES and append to the list
                smiles = Chem.MolToSmiles(mol)
                file_names.append(file_name_without_extension)
                smiles_list.append(smiles)

    return file_names, smiles_list

def get_sequence(file):
    sequence = ''
    with open(file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                continue  # Skip the header line
            sequence += line.strip()  # Append the sequence removing any whitespace
    return sequence