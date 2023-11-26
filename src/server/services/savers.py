import os
import json
from pathlib import Path

from werkzeug.datastructures import FileStorage


class FileSaver:
    def save(self, content, folder, filename):
        raise NotImplementedError

class PDBFileSaver(FileSaver):
    def save(self, content, folder, filename):
        return self._save_content(content, folder, filename)

    def pdb_to_fasta(self, pdb_filename, pdb_content, fasta_filename):
        """
        Convert a PDB file to a FASTA file.
        
        :param pdb_filename: The path to the input PDB file.
        :param fasta_filename: The path to the output FASTA file.
        """

        sequence = ""
        current_chain = None

        for line in pdb_content:
            if line.startswith("ATOM") and line[13:15].strip() == "CA":
                chain_id = line[21]
                amino_acid = line[17:20].strip()

                # If a new chain starts, separate it with a newline in the FASTA file
                if current_chain and current_chain != chain_id:
                    sequence += "\n"
                current_chain = chain_id

                # Convert three-letter amino acid codes to single-letter codes
                amino_acid_code = self.three_to_one(amino_acid)
                sequence += amino_acid_code

        with open(fasta_filename, 'w') as fasta_file:
            fasta_file.write(">Converted from {}\n".format(pdb_filename))
            fasta_file.write(sequence)

        return fasta_filename


    def three_to_one(self, three_letter_code):
        """
        Convert a three-letter amino acid code to a single-letter code.
        """
        conversion = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
            "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G",
            "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
            "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
            "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
        }
        return conversion.get(three_letter_code, "?")

    def _save_content(self, content, folder, filename):
        if not filename.endswith('.pdb'):
            filename += '.pdb'
        if not os.path.exists(folder):
            os.makedirs(folder)
        path = os.path.join(folder, filename)
        if isinstance(content, FileStorage):
            content.save(path)
            return
        with open(path, 'w') as f:
            f.write(content)

class CSVFileSaver(FileSaver):
    def save(self, content, folder, filename):
        self._save_content(content, folder, filename)

    def _save_content(self, content, folder, filename):
        if not filename.endswith('.csv'):
            filename += '.csv'
        if not os.path.exists(folder):
            os.makedirs(folder)
        with open(os.path.join(folder, filename), 'w') as f:
            f.write(content)

class SDFFileSaver(FileSaver):
    def save(self, content, folder, filename):
        return self._save_content(content, folder, filename)

    def _save_content(self, content, folder, filename):
        if not filename.endswith('.sdf'):
            filename += '.sdf'
        if not os.path.exists(folder):
            os.makedirs(folder)
        path = os.path.join(folder, filename)
        if isinstance(content, FileStorage):
            content.save(path)
            return path
        with open(path, 'w') as f:
            f.write(content)
        return path

class JSONFileSaver(FileSaver):
    def save(self, content, folder, filename):
        if not isinstance(content, dict):
            raise ValueError("The content must be a dictionary for JSONFileSaver.")
        self._save_content(content, folder, filename)

    def _save_content(self, content, folder, filename):
        if not os.path.exists(folder):
            os.makedirs(folder)
        if not filename.endswith('.json'):
            filename += '.json'
        with open(os.path.join(folder, filename), 'w') as f:
            json.dump(content, f, indent=4, default=float)

class FastaFileSaver(FileSaver):
    def save(self, content, folder, filename):
        return self._save_content(content, folder, filename)

    def _save_content(self, content, folder, filename):
        if not filename.endswith('.fasta'):
            filename += '.fasta'
        if not os.path.exists(folder):
            os.makedirs(folder)
        path = os.path.join(folder, filename)
        if isinstance(content, FileStorage):
            content.save(path)
            return path
        with open(path, 'w') as f:
            f.write(content)
        return path

class FileSaverFactory:
    @staticmethod
    def get_saver(filename):
        file_extension = os.path.splitext(filename)[1].lower()

        if file_extension == '.pdb':
            return PDBFileSaver()
        elif file_extension == '.csv':
            return CSVFileSaver()
        elif file_extension == '.sdf':
            return SDFFileSaver()
        elif file_extension == '.json':
            return JSONFileSaver()
        elif file_extension == 'skip':
            pass
        else:
            raise ValueError(f"Unsupported file extension: {file_extension}")
