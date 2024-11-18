import os
import tempfile

from rdkit import Chem
from rdkit.Chem import AllChem


def smiles_to_sdf_string(input_smiles: str) -> str:
    temp_file_descriptor, temp_file_path = tempfile.mkstemp()
    try:
        with Chem.SDWriter(temp_file_path) as sdf_writer:
            molecule = Chem.MolFromSmiles(input_smiles)
            if molecule is not None:
                AllChem.Compute2DCoords(molecule)
                sdf_writer.write(molecule)
        with open(temp_file_path, "r") as file:
            sdf_contents = file.read()
    finally:
        # Ensure the temporary file is removed
        os.close(temp_file_descriptor)
        os.remove(temp_file_path)

    return sdf_contents
