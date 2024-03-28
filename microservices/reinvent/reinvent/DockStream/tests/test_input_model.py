import unittest
from glob import glob

from dockstream.core.input_model import AzdockInput
from dockstream.utils.files_paths import attach_root_path


class TestAllJsonFilesInExamples(unittest.TestCase):
    def test_docking_jsons(self):
        files_dir = attach_root_path("examples/docking")
        json_files = glob(f"{files_dir}/*.json")
        for f in json_files:
            with self.subTest(f=f):
                inp = AzdockInput.parse_file(f)
                self.assertTrue(isinstance(inp, AzdockInput))

    def test_ligand_jsons(self):
        files_dir = attach_root_path("examples/ligand_preparation")
        json_files = glob(f"{files_dir}/*.json")
        for f in json_files:
            with self.subTest(f=f):
                inp = AzdockInput.parse_file(f)
                self.assertTrue(isinstance(inp, AzdockInput))


if __name__ == '__main__':
    unittest.main()
