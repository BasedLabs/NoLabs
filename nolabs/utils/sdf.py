class SDFReader:
    def __init__(self):
        pass

    def read_sdf(self, file_path: str) -> str:
        """
        Read and return the content of a SDF file as a string.
        """
        with open(file_path, 'r') as file:
            return file.read()


class SDFWriter:
    def __init__(self):
        pass

    def write_sdf(self, sdf_content: str, file_path: str) -> None:
        """
        Write the given SDF content to a file.
        """
        with open(file_path, 'w') as file:
            file.write(sdf_content)
