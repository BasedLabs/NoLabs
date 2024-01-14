from typing import List

from rdkit import Chem

class SDFReader:
    def __init__(self):
        pass

    def read_sdf(self, file_path: str) -> str:
        """
        Read and return the content of an SDF file as a string.
        """
        with open(file_path, 'r') as file:
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
        with open(file_path, 'w') as file:
            file.write(sdf_content)
