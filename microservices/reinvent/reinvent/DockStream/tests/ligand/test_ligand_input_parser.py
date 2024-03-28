import unittest

from dockstream.core.ligand.ligand_input_parser import LigandInputParser
from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum
from dockstream.utils.enums.docking_enum import DockingConfigurationEnum

from tests.tests_paths import PATHS_1UYD
from dockstream.utils.files_paths import attach_root_path

_LP = LigandPreparationEnum()
_CE = DockingConfigurationEnum()


class Test_ligand_input_parser(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

    def test_load_smi(self):
        conf = {_LP.INPUT: {_LP.INPUT_PATH: attach_root_path(PATHS_1UYD.LIGANDS_SMILES_TXT),
                            _LP.INPUT_TYPE: _LP.INPUT_TYPE_SMI}}
        parser = LigandInputParser(**conf)
        self.assertEqual(len(parser.get_ligands()), 15)
        self.assertEqual("C#CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2Cl)nc2c(N)ncnc21", parser.get_ligands()[0].get_smile())

        # inferring type
        conf = {_LP.INPUT: {_LP.INPUT_PATH: attach_root_path(PATHS_1UYD.LIGANDS_SMILES_SMI_TAUT_ENUM)}}
        parser = LigandInputParser(**conf)
        self.assertEqual(len(parser.get_ligands()), 17)
        self.assertEqual("C#CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2Cl)nc2c(N)ncnc21", parser.get_ligands()[0].get_smile())

    def test_load_sdf(self):
        # without tags
        conf = {_LP.INPUT: {_LP.INPUT_PATH: attach_root_path(PATHS_1UYD.LIGANDS_WITH_ENUMERATION_SDF),
                            _LP.INPUT_TYPE: _LP.INPUT_TYPE_SDF}}
        parser = LigandInputParser(**conf)

        self.assertEqual(len(parser.get_ligands()), 21)
        self.assertEqual("[H]C#CC([H])([H])C([H])([H])C([H])([H])n1c(C([H])([H])c2c([H])c(OC([H])([H])[H])c(OC([H])([H])[H])c(OC([H])([H])[H])c2Cl)nc2c(N([H])[H])nc([H])nc21",
                         parser.get_ligands()[0].get_smile())
        self.assertEqual("20:0",
                         parser.get_ligands()[20].get_identifier())

        # without tags - inferring type
        conf = {_LP.INPUT: {_LP.INPUT_PATH: attach_root_path(PATHS_1UYD.LIGANDS_WITH_ENUMERATION_SDF)}}
        parser = LigandInputParser(**conf)

        self.assertEqual(len(parser.get_ligands()), 21)
        self.assertEqual("[H]C#CC([H])([H])C([H])([H])C([H])([H])n1c(C([H])([H])c2c([H])c(OC([H])([H])[H])c(OC([H])([H])[H])c(OC([H])([H])[H])c2Cl)nc2c(N([H])[H])nc([H])nc21",
                         parser.get_ligands()[0].get_smile())

        # with tags
        conf = {_LP.INPUT: {_LP.INPUT_PATH: attach_root_path(PATHS_1UYD.LIGANDS_SDF),
                            _LP.INPUT_TYPE: _LP.INPUT_TYPE_SDF,
                            _LP.INPUT_SDF_TAGS: {_LP.INPUT_SDF_TAGNAME_NAMES: "name"}}}
        parser = LigandInputParser(**conf)

        self.assertEqual(len(parser.get_ligands()), 15)
        self.assertEqual("C#CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2Cl)nc2c(N)ncnc21",
                         parser.get_ligands()[0].get_smile())
        self.assertEqual(parser.get_ligands()[0].get_name(), "1uye-A-1224 PU9 Conf:1")

        # without tags, initialization mode "dockstream" -> groups enumerations according to DockStream nomenclature
        conf = {_LP.INPUT: {_LP.INPUT_PATH: attach_root_path(PATHS_1UYD.LIGANDS_WITH_ENUMERATION_SDF),
                            _LP.INITIALIZATION_MODE: _LP.INITIALIZATION_MODE_AZDOCK,
                            _LP.INPUT_TYPE: _LP.INPUT_TYPE_SDF}}
        parser = LigandInputParser(**conf)

        self.assertEqual(len(parser.get_ligands()), 21)
        self.assertEqual("[H]C#CC([H])([H])C([H])([H])C([H])([H])n1c(C([H])([H])c2c([H])c(OC([H])([H])[H])c(OC([H])([H])[H])c(OC([H])([H])[H])c2Cl)nc2c(N([H])[H])nc([H])nc21",
                         parser.get_ligands()[0].get_smile())
        self.assertEqual("16:3",
                         parser.get_ligands()[20].get_identifier())

    def test_load_csv(self):
        # without "names"
        conf = {_LP.INPUT: {_LP.INPUT_PATH: attach_root_path(PATHS_1UYD.LIGANDS_CSV),
                            _LP.INPUT_TYPE: _LP.INPUT_TYPE_CSV,
                            _LP.INPUT_CSV_COLUMNS: {
                                _LP.INPUT_CSV_COLNAME_SMILES: "smiles"
                            }}}
        parser = LigandInputParser(**conf)

        self.assertEqual(len(parser.get_ligands()), 17)
        self.assertEqual("C#CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2Cl)nc2c(N)ncnc21", parser.get_ligands()[0].get_smile())
        self.assertIsNone(parser.get_ligands()[16].get_name())

        # without "names" - inferring type
        conf = {_LP.INPUT: {_LP.INPUT_PATH: attach_root_path(PATHS_1UYD.LIGANDS_CSV),
                            _LP.INPUT_CSV_COLUMNS: {
                                _LP.INPUT_CSV_COLNAME_SMILES: "smiles"}}}
        parser = LigandInputParser(**conf)

        self.assertEqual(len(parser.get_ligands()), 17)
        self.assertEqual("C#CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2Cl)nc2c(N)ncnc21", parser.get_ligands()[0].get_smile())
        self.assertIsNone(parser.get_ligands()[16].get_name())

        # with "names"
        conf = {_LP.INPUT: {_LP.INPUT_PATH: attach_root_path(PATHS_1UYD.LIGANDS_CSV),
                            _LP.INPUT_TYPE: _LP.INPUT_TYPE_CSV,
                            _LP.INPUT_CSV_COLUMNS: {
                                _LP.INPUT_CSV_COLNAME_SMILES: "smiles",
                                _LP.INPUT_CSV_COLNAME_NAMES: "name"
                }}}
        parser = LigandInputParser(**conf)

        self.assertEqual(len(parser.get_ligands()), 17)
        self.assertEqual("C#CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2Cl)nc2c(N)ncnc21", parser.get_ligands()[0].get_smile())
        self.assertEqual(parser.get_ligands()[16].get_name(), "taut_enum_comp2")

    def test_load_console(self):
        conf = {_LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_CONSOLE}}
        parser = LigandInputParser(smiles="C#CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2Cl)nc2c(N)ncnc21;CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2)nc2c(N)ncnc21;CCCCn1c(Cc2cc(OC)ccc2OC)nc2c(N)ncnc21",
                                   **conf)

        self.assertEqual(len(parser.get_ligands()), 3)
        self.assertEqual("C#CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2Cl)nc2c(N)ncnc21", parser.get_ligands()[0].get_smile())

        # inferring type
        conf = {_LP.INPUT: {}}
        parser = LigandInputParser(
            smiles="C#CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2Cl)nc2c(N)ncnc21;CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2)nc2c(N)ncnc21;CCCCn1c(Cc2cc(OC)ccc2OC)nc2c(N)ncnc21",
            **conf)

        self.assertEqual(len(parser.get_ligands()), 3)
        self.assertEqual("C#CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2Cl)nc2c(N)ncnc21", parser.get_ligands()[0].get_smile())
