import unittest

from dockstream.core.TautEnum.taut_enum_smile_preparation import TautEnumSmilePreparator

from dockstream.utils.enums.docking_enum import DockingConfigurationEnum
from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum
from dockstream.utils.enums.taut_enum_enums import TautEnumEnum

from tests.tests_paths import PATHS_1UYD
from dockstream.utils.files_paths import attach_root_path
from dockstream.utils.smiles import read_smiles_file
from dockstream.core.ligand.ligand import Ligand


# TODO: move "prefix_execution" to "config/tests_config/config.json"
class Test_Corina_ligand_preparation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._DE = DockingConfigurationEnum()
        cls._LP = LigandPreparationEnum()
        cls._TE = TautEnumEnum()

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

    @classmethod
    def tearDownClass(cls):
        pass

    def test_taut_enum_with_protonation(self):
        prep = TautEnumSmilePreparator(enumerate_protonation=True,
                                       original_enumeration=True,
                                       add_numbers_to_name=True,
                                       prefix_execution="module load taut_enum")
        taut_result = prep.annotate_tautomers(ligands=self.ligands)

        # check tautomer assignment
        list_smiles = []
        list_ligand_numbers = []
        list_enumeration_numbers = []
        for result in taut_result:
            list_smiles.append(result.get_smile())
            list_ligand_numbers.append(result.get_ligand_number())
            list_enumeration_numbers.append(result.get_enumeration())
        self.assertListEqual(list_smiles,
                             ['COc1cc(c(c(c1OC)OC)Cl)Cc2nc3c(ncnc3n2CCCC#C)N',
                              'CCCCn1c(nc2c1ncnc2N)Cc3cc(c(c(c3)OC)OC)OC', 'CCCCn1c(nc2c1ncnc2N)Cc3cc(ccc3OC)OC',
                              'CCCCn1c(nc2c1ncnc2N)Cc3cccc(c3)OC', 'COc1cc(c(c(c1OC)OC)Cl)Cc2nc3c(nc(nc3n2CCCC#C)F)N',
                              'CCCCn1c(nc2c1ncnc2N)Cc3ccc(cc3)OC', 'CCCCn1c(nc2c1ncnc2N)Cc3ccc4c(c3)OCO4',
                              'CCCCn1c(nc2c1nc(nc2N)F)Cc3cc(ccc3OC)OC', 'CCCCn1c(nc2c1nc(nc2N)F)Cc3ccc4c(c3)OCO4',
                              'COc1ccc(c(c1)Cc2nc3c(nc(nc3n2CCCC#C)F)N)OC',
                              'CC(C)[NH2+]CCCn1c(nc2c1nc(nc2N)F)Cc3cc4c(cc3I)OCO4',
                              'CC(C)[NH2+]CCCn1c2c(c(ncn2)N)nc1Sc3cc4c(cc3Br)OCO4',
                              'CC(C)[NH2+]CCCn1c2c(c(ncn2)N)nc1Sc3cc4c(cc3I)OCO4',
                              'COc1ccc(c(c1)Cc2[nH]c3c(n2)c(nc(n3)F)N)OC',
                              'c1ccc(cc1)CNc2c(nc3n2ccnc3N)Cc4cc5c(cc4Br)OCO5', 'Cc1cc(cn1C)c2csc(n2)[NH+]=C(N)N',
                              'Cc1cc(cn1C)c2csc(n2)N=C(N)N',
                              'CC(C)(C)CCC1(c2ccccc2C(=O)[C-](C1=O)C3=NS(=O)(=O)c4cc(ccc4[N-]3)[N-]S(=O)(=O)C)[NH2+]Cc5cccc(c5)OC',
                              'CC(C)(C)CCC1(c2ccccc2C(=O)[C-](C1=O)C3=NS(=O)(=O)c4cc(ccc4[N-]3)NS(=O)(=O)C)[NH2+]Cc5cccc(c5)OC',
                              'CC(C)(C)CCC1(c2ccccc2C(=O)[C-](C1=O)C3=NS(=O)(=O)c4cc(ccc4N3)[N-]S(=O)(=O)C)[NH2+]Cc5cccc(c5)OC',
                              'CC(C)(C)CCC1(c2ccccc2C(=O)[C-](C1=O)C3=NS(=O)(=O)c4cc(ccc4N3)NS(=O)(=O)C)[NH2+]Cc5cccc(c5)OC'])
        self.assertListEqual(list_ligand_numbers, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 15, 16, 16, 16, 16])
        self.assertListEqual(list_enumeration_numbers, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 2, 3])

    def test_taut_enum_without_protonation(self):
        prep = TautEnumSmilePreparator(enumerate_protonation=False,
                                       original_enumeration=True,
                                       add_numbers_to_name=True,
                                       prefix_execution="module load taut_enum",
                                       binary_location=None)
        taut_result = prep.annotate_tautomers(ligands=self.ligands)

        # check tautomer assignment
        list_smiles = []
        list_ligand_numbers = []
        list_enumeration_numbers = []
        for result in taut_result:
            list_smiles.append(result.get_smile())
            list_ligand_numbers.append(result.get_ligand_number())
            list_enumeration_numbers.append(result.get_enumeration())
        self.assertListEqual(list_smiles,
                             ['COc1cc(c(c(c1OC)OC)Cl)Cc2nc3c(ncnc3n2CCCC#C)N',
                              'CCCCn1c(nc2c1ncnc2N)Cc3cc(c(c(c3)OC)OC)OC', 'CCCCn1c(nc2c1ncnc2N)Cc3cc(ccc3OC)OC',
                              'CCCCn1c(nc2c1ncnc2N)Cc3cccc(c3)OC', 'COc1cc(c(c(c1OC)OC)Cl)Cc2nc3c(nc(nc3n2CCCC#C)F)N',
                              'CCCCn1c(nc2c1ncnc2N)Cc3ccc(cc3)OC', 'CCCCn1c(nc2c1ncnc2N)Cc3ccc4c(c3)OCO4',
                              'CCCCn1c(nc2c1nc(nc2N)F)Cc3cc(ccc3OC)OC', 'CCCCn1c(nc2c1nc(nc2N)F)Cc3ccc4c(c3)OCO4',
                              'COc1ccc(c(c1)Cc2nc3c(nc(nc3n2CCCC#C)F)N)OC',
                              'CC(C)NCCCn1c(nc2c1nc(nc2N)F)Cc3cc4c(cc3I)OCO4',
                              'CC(C)NCCCn1c2c(c(ncn2)N)nc1Sc3cc4c(cc3Br)OCO4',
                              'CC(C)NCCCn1c2c(c(ncn2)N)nc1Sc3cc4c(cc3I)OCO4',
                              'COc1ccc(c(c1)Cc2[nH]c3c(n2)c(nc(n3)F)N)OC',
                              'c1ccc(cc1)CNc2c(nc3n2ccnc3N)Cc4cc5c(cc4Br)OCO5', 'Cc1cc(cn1C)c2csc(n2)NC(=N)N',
                              'Cc1cc(cn1C)c2csc(n2)N=C(N)N',
                              'CC(C)(C)CCC1(c2ccccc2C(=O)C(C1=O)C3=Nc4ccc(cc4S(=O)(=O)N3)NS(=O)(=O)C)NCc5cccc(c5)OC',
                              'CC(C)(C)CCC1(c2ccccc2C(=O)C(C1=O)C3=NS(=O)(=O)c4cc(ccc4N3)NS(=O)(=O)C)NCc5cccc(c5)OC'])

        self.assertListEqual(list_ligand_numbers, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 15, 16, 16])
        self.assertListEqual(list_enumeration_numbers, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1])