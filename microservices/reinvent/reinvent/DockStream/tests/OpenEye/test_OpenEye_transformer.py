import unittest
import os

from dockstream.utils.enums.transformations_enums import TransformationEnum
from dockstream.utils.enums.OpenEye_enums import OpenEyeDockingConfigurationEnum, OpenEyeLigandPreparationEnum

from dockstream.core.ligand.ligand import Ligand
from dockstream.core.OpenEye.OpenEye_transformator import OpenEyeTransformator

from tests.tests_paths import PATHS_1UYD, MAIN_CONFIG
from dockstream.utils.files_paths import attach_root_path
from dockstream.utils.smiles import read_smiles_file


class Test_OpenEye_transformer(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._LP = OpenEyeLigandPreparationEnum()
        cls._CE = OpenEyeDockingConfigurationEnum()
        cls._TE = TransformationEnum()
        if "OE_LICENSE" in MAIN_CONFIG:
            os.environ["OE_LICENSE"] = MAIN_CONFIG["OE_LICENSE"]

    def setUp(self):
        self.ligands_smiles = list(read_smiles_file(attach_root_path(PATHS_1UYD.LIGANDS_SMILES_TXT), standardize=False))
        ligands = []
        for smile_number, smile in enumerate(self.ligands_smiles):
            ligands.append(Ligand(smile=smile, ligand_number=smile_number))
        self.ligands = ligands

    @classmethod
    def tearDownClass(cls):
        pass

    def test_smirk_transformation_keep(self):
        # use a smirk, that will attach a charge and a hydrogen atom to a given nitrogen
        conf = {self._TE.TRANSFORMATION_TYPE: self._TE.TRANSFORMATION_TYPE_SMIRKS,
                self._TE.TRANSFORMATION_BACKEND: self._TE.TRANSFORMATION_BACKEND_OPENEYE,
                self._TE.TRANSFORMATION_SMIRKS: "[c:1]1[n:2][c:3]2[c:4]([n:5][c:6][n;X2:7][c:8]2[n:9]1)>>[c:1]1[n:2][c:3]2[c:4]([n:5][c:6][n+:7]([H])[c:8]2[n:9]1)",
                self._TE.TRANSFORMATION_FAIL_ACTION: self._TE.TRANSFORMATION_FAIL_ACTION_KEEP}

        OE_transformator = OpenEyeTransformator(conf)
        self.assertListEqual(["C#CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2Cl)nc2c(N)ncnc21",
                              "COc1ccc(OC)c(Cc2nc3nc(F)nc(N)c3[nH]2)c1",
                              "Nc1nccn2c(NCc3ccccc3)c(Cc3cc4c(cc3Br)OCO4)nc12"],
                             [self.ligands[index].get_smile() for index in [0, 13, 14]])

        # the last smile does not have match the substructure and as "fail_action" is "keep" will remain untouched
        transformed_ligands = OE_transformator.transform(self.ligands)
        self.assertEqual(len(transformed_ligands), 15)
        self.assertListEqual(["COc1cc(c(c(c1OC)OC)Cl)Cc2nc3c(nc[nH+]c3n2CCCC#C)N",
                              "COc1ccc(c(c1)Cc2[nH]c3c(nc([nH+]c3n2)F)N)OC",
                              "Nc1nccn2c(NCc3ccccc3)c(Cc3cc4c(cc3Br)OCO4)nc12"],
                             [self.ligands[index].get_smile() for index in [0, 13, 14]])

    def test_smirk_transformation_discard(self):
        # use a smirk, that will attach a charge and a hydrogen atom to a given nitrogen
        conf = {self._TE.TRANSFORMATION_TYPE: self._TE.TRANSFORMATION_TYPE_SMIRKS,
                self._TE.TRANSFORMATION_BACKEND: self._TE.TRANSFORMATION_BACKEND_OPENEYE,
                self._TE.TRANSFORMATION_SMIRKS: "[c:1]1[n:2][c:3]2[c:4]([n:5][c:6][n;X2:7][c:8]2[n:9]1)>>[c:1]1[n:2][c:3]2[c:4]([n:5][c:6][n+:7]([H])[c:8]2[n:9]1)",
                self._TE.TRANSFORMATION_FAIL_ACTION: self._TE.TRANSFORMATION_FAIL_ACTION_DISCARD}

        OE_transformator = OpenEyeTransformator(conf)
        self.assertListEqual(["C#CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2Cl)nc2c(N)ncnc21",
                              "COc1ccc(OC)c(Cc2nc3nc(F)nc(N)c3[nH]2)c1",
                              "Nc1nccn2c(NCc3ccccc3)c(Cc3cc4c(cc3Br)OCO4)nc12"],
                             [self.ligands[index].get_smile() for index in [0, 13, 14]])

        # the last smile does not have match the substructure and as "fail_action" is "keep" will remain untouched
        transformed_ligands = OE_transformator.transform(self.ligands)
        self.assertEqual(len(transformed_ligands), 14)
        self.assertListEqual(["COc1cc(c(c(c1OC)OC)Cl)Cc2nc3c(nc[nH+]c3n2CCCC#C)N",
                              "COc1ccc(c(c1)Cc2[nH]c3c(nc([nH+]c3n2)F)N)OC"],
                             [self.ligands[index].get_smile() for index in [0, 13]])
