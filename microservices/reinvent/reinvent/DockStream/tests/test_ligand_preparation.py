import unittest

from rdkit import Chem

from dockstream.core.ligand_preparator import LigandPreparator, Input
from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum
from dockstream.utils.enums.docking_enum import DockingConfigurationEnum

from tests.tests_paths import PATHS_1UYD
from dockstream.utils.files_paths import attach_root_path

_LP = LigandPreparationEnum()
_CE = DockingConfigurationEnum()


class Test_ligand_preparation(unittest.TestCase):

    def setUp(self):
        file = attach_root_path(PATHS_1UYD.LIGANDS_SDF)
        mols = [mol for mol in Chem.SDMolSupplier(file) if mol is not None]
        self.ligands = mols

    def test_initialization(self):
        prep = LigandPreparator(
            ligands=self.ligands,
            pool_id="testPool",
            input=Input()
        )
        self.assertEqual(len(prep.get_ligands()), 15)
        self.assertEqual(type(prep.get_ligands()), list)
