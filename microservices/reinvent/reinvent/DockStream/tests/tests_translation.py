import os
import unittest

from rdkit import Chem
import openeye.oechem as oechem

from tests.tests_paths import PATHS_1UYD, MAIN_CONFIG
from dockstream.utils.files_paths import attach_root_path
from dockstream.utils.translations.translation import OpenEyeMolToRDkitMol, RDkitMolToOpenEyeMol
from dockstream.utils.smiles import to_smiles


class Test_molecule_container_translation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if "OE_LICENSE" in MAIN_CONFIG:
            os.environ["OE_LICENSE"] = MAIN_CONFIG["OE_LICENSE"]

    def setUp(self):
        # read in the molecules in "OpenEye" containers
        ifs = oechem.oemolistream()
        ifs.SetFormat(oechem.OEFormat_SDF)
        ifs.open(attach_root_path(PATHS_1UYD.LIGANDS_SDF))

        self.ligands_OpenEye = []
        for mol in ifs.GetOEMols():
            # this is critical: make a copy, otherwise empty molecules will populate the list
            self.ligands_OpenEye.append(oechem.OEMol(mol))

        # read in the molecules in "RDkit" format
        self.ligands_RDkit = []
        for mol in Chem.SDMolSupplier(attach_root_path(PATHS_1UYD.LIGANDS_SDF)):
            if mol is None:
                continue
            self.ligands_RDkit.append(mol)

    @classmethod
    def tearDownClass(cls):
        pass

    def test_RDkit2OpenEye(self):
        translated_ligands_OpenEye = []
        for mol in self.ligands_RDkit:
            translated_ligands_OpenEye.append(RDkitMolToOpenEyeMol(molecule=mol, bySMILES=False))
        self.assertEqual(len(translated_ligands_OpenEye), 15)

        m1_trans = translated_ligands_OpenEye[1]
        m1_ori = self.ligands_OpenEye[1]
        self.assertEqual([m1_ori.NumAtoms(),
                          oechem.OECount(m1_ori, oechem.OEIsHeavy()),
                          oechem.OECount(m1_ori, oechem.OEAtomIsInRing()),
                          oechem.OECount(m1_ori, oechem.OEIsRotor())],
                         [m1_trans.NumAtoms(),
                          oechem.OECount(m1_trans, oechem.OEIsHeavy()),
                          oechem.OECount(m1_trans, oechem.OEAtomIsInRing()),
                          oechem.OECount(m1_trans, oechem.OEIsRotor())])

    def test_OpenEye2RDkit(self):
        translated_ligands_RDkit = []
        for mol in self.ligands_OpenEye:
            translated_ligands_RDkit.append(OpenEyeMolToRDkitMol(molecule=mol, bySMILES=False))
        self.assertEqual(len(translated_ligands_RDkit), 15)

        m1_trans = translated_ligands_RDkit[1]
        m1_ori = self.ligands_RDkit[1]
        self.assertEqual([m1_ori.GetNumAtoms(), m1_ori.GetNumHeavyAtoms(), m1_ori.GetNumBonds()],
                         [m1_trans.GetNumAtoms(), m1_trans.GetNumHeavyAtoms(), m1_trans.GetNumBonds()])
        self.assertEqual(to_smiles(m1_ori),
                         to_smiles(m1_trans))

    def test_RDkit2OpenEye_bySMILES(self):
        translated_ligands_OpenEye = []
        for mol in self.ligands_RDkit:
            translated_ligands_OpenEye.append(RDkitMolToOpenEyeMol(molecule=mol, bySMILES=True))
        self.assertEqual(len(translated_ligands_OpenEye), 15)

        m1_trans = translated_ligands_OpenEye[1]
        m1_ori = self.ligands_OpenEye[1]
        self.assertEqual([m1_ori.NumAtoms(),
                          oechem.OECount(m1_ori, oechem.OEIsHeavy()),
                          oechem.OECount(m1_ori, oechem.OEAtomIsInRing()),
                          oechem.OECount(m1_ori, oechem.OEIsRotor())],
                         [m1_trans.NumAtoms(),
                          oechem.OECount(m1_trans, oechem.OEIsHeavy()),
                          oechem.OECount(m1_trans, oechem.OEAtomIsInRing()),
                          oechem.OECount(m1_trans, oechem.OEIsRotor())])
        self.assertEqual(oechem.OEMolToSmiles(m1_ori),
                         oechem.OEMolToSmiles(m1_trans))

    def test_OpenEye2RDkit_bySMILES(self):
        translated_ligands_RDkit = []
        for mol in self.ligands_OpenEye:
            translated_ligands_RDkit.append(OpenEyeMolToRDkitMol(molecule=mol, bySMILES=True))
        self.assertEqual(len(translated_ligands_RDkit), 15)

        m1_trans = translated_ligands_RDkit[1]
        m1_ori = self.ligands_RDkit[1]
        self.assertEqual([m1_ori.GetNumAtoms(), m1_ori.GetNumHeavyAtoms(), m1_ori.GetNumBonds()],
                         [m1_trans.GetNumAtoms(), m1_trans.GetNumHeavyAtoms(), m1_trans.GetNumBonds()])
        self.assertEqual(to_smiles(m1_ori),
                         to_smiles(m1_trans))
