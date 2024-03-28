import unittest

from dockstream.core.RDkit.RDkit_stereo_enumerator import RDKitStereoEnumerator, RDKitStereoEnumeratorParameters

from tests.tests_paths import PATHS_1UYD
from dockstream.utils.files_paths import attach_root_path
from dockstream.utils.smiles import read_smiles_file
from dockstream.core.ligand.ligand import Ligand


# TODO: move "prefix_execution" to "config/tests_config/config.json"
class Test_RDkit_stereo_enumeration(unittest.TestCase):

        def setUp(self):
            smiles = list(read_smiles_file(attach_root_path(PATHS_1UYD.LIGANDS_SMILES_TXT), standardize=False))

            # add one smiles that will result in more than one enumeration (unlike the 15 others)
            smiles.append("Cc1cc(cn1C)c2csc(n2)N=C(N)N")
            smiles.append("CC(C)(C)CCC1(c2ccccc2C(=O)C(C1=O)C3=NS(=O)(=O)c4cc(ccc4N3)NS(=O)(=O)C)NCc5cccc(c5)OC")

            self.ligands = []
            for smile_number, smile in enumerate(smiles):
                self.ligands.append(Ligand(smile=smile,
                                           original_smile=smile,
                                           ligand_number=smile_number,
                                           enumeration=0,
                                           molecule=None,
                                           mol_type=None))

        def test_RDkit_stereo_enumerator(self):
            enum = RDKitStereoEnumerator(
                parameters=RDKitStereoEnumeratorParameters(
                    try_embedding=True,
                    unique=True,
                    max_isomers=1024))
            enum_result = enum.enumerate(ligands=self.ligands)

            # check tautomer assignment
            list_smiles = []
            list_ligand_numbers = []
            list_enumeration_numbers = []
            for result in enum_result:
                list_smiles.append(result.get_smile())
                list_ligand_numbers.append(result.get_ligand_number())
                list_enumeration_numbers.append(result.get_enumeration())
            self.assertListEqual(list_smiles,
                                 ['C#CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2Cl)nc2c(N)ncnc21',
                                  'CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2)nc2c(N)ncnc21',
                                  'CCCCn1c(Cc2cc(OC)ccc2OC)nc2c(N)ncnc21',
                                  'CCCCn1c(Cc2cccc(OC)c2)nc2c(N)ncnc21',
                                  'C#CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2Cl)nc2c(N)nc(F)nc21',
                                  'CCCCn1c(Cc2ccc(OC)cc2)nc2c(N)ncnc21',
                                  'CCCCn1c(Cc2ccc3c(c2)OCO3)nc2c(N)ncnc21',
                                  'CCCCn1c(Cc2cc(OC)ccc2OC)nc2c(N)nc(F)nc21',
                                  'CCCCn1c(Cc2ccc3c(c2)OCO3)nc2c(N)nc(F)nc21',
                                  'C#CCCCn1c(Cc2cc(OC)ccc2OC)nc2c(N)nc(F)nc21',
                                  'CC(C)NCCCn1c(Cc2cc3c(cc2I)OCO3)nc2c(N)nc(F)nc21',
                                  'CC(C)NCCCn1c(Sc2cc3c(cc2Br)OCO3)nc2c(N)ncnc21',
                                  'CC(C)NCCCn1c(Sc2cc3c(cc2I)OCO3)nc2c(N)ncnc21',
                                  'COc1ccc(OC)c(Cc2nc3nc(F)nc(N)c3[nH]2)c1',
                                  'Nc1nccn2c(NCc3ccccc3)c(Cc3cc4c(cc3Br)OCO4)nc12',
                                  'Cc1cc(-c2csc(N=C(N)N)n2)cn1C',
                                  'COc1cccc(CN[C@@]2(CCC(C)(C)C)C(=O)[C@@H](C3=NS(=O)(=O)c4cc(NS(C)(=O)=O)ccc4N3)C(=O)c3ccccc32)c1',
                                  'COc1cccc(CN[C@@]2(CCC(C)(C)C)C(=O)[C@H](C3=NS(=O)(=O)c4cc(NS(C)(=O)=O)ccc4N3)C(=O)c3ccccc32)c1',
                                  'COc1cccc(CN[C@]2(CCC(C)(C)C)C(=O)[C@@H](C3=NS(=O)(=O)c4cc(NS(C)(=O)=O)ccc4N3)C(=O)c3ccccc32)c1',
                                  'COc1cccc(CN[C@]2(CCC(C)(C)C)C(=O)[C@H](C3=NS(=O)(=O)c4cc(NS(C)(=O)=O)ccc4N3)C(=O)c3ccccc32)c1'])
            self.assertListEqual(list_ligand_numbers,
                                 [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 16, 16, 16])
            self.assertListEqual(list_enumeration_numbers,
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3])
