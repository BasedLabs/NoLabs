import os
import tempfile

from rdkit import Chem
from rdkit.Chem import AllChem


class SDFReader:
    def __init__(self):
        pass

    def read_sdf(self, file_path: str) -> str:
        """
        Read and return the content of an SDF file as a string.
        """
        with open(file_path, "r") as file:
            return file.read()

    def get_smiles_from_sdf(self, file_path: str) -> str:
        """
        Read and return the smiles string of an SDF file. (now it reads only the first molecule)
        """
        sdf_data = Chem.SDMolSupplier(file_path)

        smiles_list = []
        for mol in sdf_data:
            if mol is not None:
                smiles = Chem.MolToSmiles(mol)
                smiles_list.append(smiles)

        return smiles_list[0]


class SDFWriter:
    def __init__(self):
        pass

    def write_sdf(self, sdf_content: str, file_path: str) -> None:
        """
        Write the given SDF content to a file.
        """
        with open(file_path, "w") as file:
            file.write(sdf_content)


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
