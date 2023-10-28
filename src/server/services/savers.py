import os
import json
from pathlib import Path

from werkzeug.datastructures import FileStorage


class FileSaver:
    def save(self, content, folder, filename):
        raise NotImplementedError

class PDBFileSaver(FileSaver):
    def save(self, content, folder, filename):
        self._save_content(content, folder, filename)

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
        self._save_content(content, folder, filename)

    def _save_content(self, content, folder, filename):
        if not filename.endswith('.sdf'):
            filename += '.sdf'
        if not os.path.exists(folder):
            os.makedirs(folder)
        path = os.path.join(folder, filename)
        if isinstance(content, FileStorage):
            content.save(path)
            return
        with open(path, 'w') as f:
            f.write(content)

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
