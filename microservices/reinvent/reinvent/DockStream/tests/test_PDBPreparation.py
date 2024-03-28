import os
import unittest
import tempfile

from rdkit import Chem

from dockstream.core.pdb_preparator import PDBPreparator
from dockstream.containers.target_preparation_container import TargetPreparationContainer
from dockstream.utils.enums.target_preparation_enum import TargetPreparationEnum
from dockstream.utils.general_utils import gen_temp_file

from tests.tests_paths import PATHS_1UYD
from dockstream.utils.files_paths import attach_root_path


class Test_PDBPreparation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._TE = TargetPreparationEnum()

    def setUp(self):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_fix_good_structure(self):
        conf = TargetPreparationContainer({self._TE.TARGETPREP: {self._TE.FIX: {self._TE.FIX_STANDARDIZE: True,
                                                                                self._TE.FIX_MISSINGHEAVYATOMS: True,
                                                                                self._TE.FIX_MISSINGLOOPS: True,
                                                                                self._TE.FIX_ADDWATERBOX: True,
                                                                                self._TE.FIX_REMOVEHETEROGENS: True,
                                                                                self._TE.FIX_MISSINGHYDROGENS: False}}})
        input_pdb_file = attach_root_path(PATHS_1UYD.TARGET_APO_PDB)
        temp_pdb_output = gen_temp_file(suffix=".pdb")
        input_mol = Chem.MolFromPDBFile(input_pdb_file, sanitize=True)

        # first processing: add water box and do not fix hydrogens that are missing
        prep = PDBPreparator(conf=conf)
        self.assertEqual(input_mol.GetNumAtoms(), 1921)
        self.assertEqual(input_mol.GetNumHeavyAtoms(), 1921)
        self.assertEqual([input_mol.GetConformer().GetPositions()[1][1],
                          input_mol.GetConformer().GetPositions()[7][1],
                          input_mol.GetConformer().GetPositions()[35][1]],
                         [28.546, 26.816, 21.825])
        prep.fix_pdb(input_pdb_file=input_pdb_file, output_pdb_file=temp_pdb_output)
        output_mol = Chem.MolFromPDBFile(temp_pdb_output, sanitize=True)
        self.assertEqual(output_mol.GetNumAtoms(), 3835)
        self.assertEqual(output_mol.GetNumHeavyAtoms(), 3835)
        os.remove(temp_pdb_output)

        # second processing: no water box and but add hydrogens
        temp_pdb_output = gen_temp_file(suffix=".pdb")
        conf[self._TE.TARGETPREP][self._TE.FIX][self._TE.FIX_ADDWATERBOX] = False
        conf[self._TE.TARGETPREP][self._TE.FIX][self._TE.FIX_MISSINGHYDROGENS] = True
        prep = PDBPreparator(conf=conf)
        prep.fix_pdb(input_pdb_file=input_pdb_file, output_pdb_file=temp_pdb_output)
        output_mol = Chem.MolFromPDBFile(temp_pdb_output, sanitize=True, removeHs=False)
        self.assertEqual(output_mol.GetNumAtoms(), 4126)
        self.assertEqual(output_mol.GetNumHeavyAtoms(), 1922)
        os.remove(temp_pdb_output)

    def test_fix_bad_structure(self):
        conf = TargetPreparationContainer({self._TE.TARGETPREP: {self._TE.FIX: {self._TE.FIX_STANDARDIZE: True,
                                                                                self._TE.FIX_MISSINGHEAVYATOMS: True,
                                                                                self._TE.FIX_MISSINGLOOPS: True,
                                                                                self._TE.FIX_ADDWATERBOX: False,
                                                                                self._TE.FIX_REMOVEHETEROGENS: True,
                                                                                self._TE.FIX_MISSINGHYDROGENS: True}}})
        input_pdb_file = attach_root_path(PATHS_1UYD.LIGAND_MISSING_PARTS_PDB)
        temp_pdb_output = gen_temp_file(suffix=".pdb")
        input_mol = Chem.MolFromPDBFile(input_pdb_file, sanitize=True)
        prep = PDBPreparator(conf=conf)
        self.assertEqual(input_mol.GetNumAtoms(), 1892)
        self.assertEqual(input_mol.GetNumHeavyAtoms(), 1892)
        prep.fix_pdb(input_pdb_file=input_pdb_file, output_pdb_file=temp_pdb_output)
        output_mol = Chem.MolFromPDBFile(temp_pdb_output, removeHs=False, proximityBonding=False, sanitize=True)
        self.assertEqual(output_mol.GetNumAtoms(), 4533)
        self.assertEqual(output_mol.GetNumHeavyAtoms(), 2144)
        os.remove(temp_pdb_output)
