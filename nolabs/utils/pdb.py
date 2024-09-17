class PDBReader:
    def __init__(self):
        pass

    def read_pdb(self, file_path: str) -> str:
        """
        Read and return the content of a PDB file as a string.
        """
        with open(file_path, "r") as file:
            return file.read()


class PDBWriter:
    def __init__(self):
        pass

    def write_pdb(self, pdb_content: str, file_path: str) -> None:
        """
        Write the given PDB content to a file.
        """
        with open(file_path, "w") as file:
            file.write(pdb_content)
