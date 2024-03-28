import unittest
import os

from dockstream.core.Corina.Corina_ligand_preparator import CorinaLigandPreparator
from dockstream.core.ligand.ligand_input_parser import LigandInputParser

from dockstream.utils.enums.docking_enum import DockingConfigurationEnum
from dockstream.utils.enums.Corina_enums import CorinaLigandPreparationEnum

from tests.tests_paths import PATHS_1UYD
from dockstream.utils.files_paths import attach_root_path
from dockstream.utils.smiles import read_smiles_file

_DE = DockingConfigurationEnum()
_LP = CorinaLigandPreparationEnum()


# TODO: move "prefix_execution" to "config/tests_config/config.json"
class Test_Corina_ligand_preparation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        self.smiles = list(read_smiles_file(attach_root_path(PATHS_1UYD.LIGANDS_SMILES_TXT),
                                            standardize=False))

    @classmethod
    def tearDownClass(cls):
        pass

    def test_coordinate_generation(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {
                    _LP.PREFIX_EXECUTION: "module load corina"}}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = CorinaLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        # check 3D coordinate generation
        prep.generate3Dcoordinates()
        self.assertEqual(prep.get_number_ligands(),
                         15)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(),
                         51)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-0.6841, -5.3454, 1.8254])

        # test write-out
        out_path = attach_root_path("tests/junk/corina_ligands.sdf")
        prep.write_ligands(path=out_path, format=_LP.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 64610)

    def test_coordinate_generation_fails(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {
                    _LP.PREFIX_EXECUTION: "module load corina"}}
        test = self.smiles
        test.append("abc")
        lig_parser = LigandInputParser(smiles=test, **conf)
        prep = CorinaLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        # check 3D coordinate generation
        self.assertEqual(prep.get_number_ligands(),
                         16)
        prep.generate3Dcoordinates()
        self.assertEqual(prep.get_number_ligands(),
                         16)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(),
                         51)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-0.6841, -5.3454, 1.8254])

        # test write-out
        out_path = attach_root_path("tests/junk/corina_ligands.sdf")
        prep.write_ligands(path=out_path, format=_LP.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 64610)

    def test_aligning(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {
                    _LP.PREFIX_EXECUTION: "module load corina"
                },
                _LP.ALIGN: {_LP.ALIGN_MODE: _LP.ALIGN_MODE_INTERNAL,
                            _LP.ALIGN_REFERENCE_PATHS: [attach_root_path(PATHS_1UYD.LIGAND_PU8_SDF)],
                            _LP.ALIGN_REFERENCE_FORMAT: _LP.ALIGN_REFERENCE_FORMAT_SDF,
                            _LP.ALIGN_MINIMUM_SUBSTRUCTURE_RATIO: 0.2,
                            _LP.ALIGN_FAIL_ACTION: _LP.ALIGN_FAIL_DISCARD,
                            _LP.ALIGN_COMPLETE_RINGS_ONLY: True}}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = CorinaLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        # check, whether alignment successfully moves the ligands
        prep.generate3Dcoordinates()
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-0.6841, -5.3454, 1.8254])
        prep.align_ligands()
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [5.274095128434708, 7.202258233600961, 21.600454578046943])

    def test_stereo_enum(self):
        smiles = self.smiles

        # add one smiles that will result in more than one enumeration (unlike the 14 others)
        smiles.append("CC(C)(C)CCC1(c2ccccc2C(=O)C(C1=O)C3=NS(=O)(=O)c4cc(ccc4N3)NS(=O)(=O)C)NCc5cccc(c5)OC")

        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {
                    _LP.PREFIX_EXECUTION: "module load corina",
                    _LP.D_OPTIONS: ["wh", "stergen", "preserve",
                                         "noflapn", "ori", "ampax",
                                         "names", "rc", "mc=1"],
                    _LP.ENUMERATE_STEREO: True
                }}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = CorinaLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        # check 3D coordinate generation
        prep.generate3Dcoordinates()
        self.assertEqual(prep.get_number_ligands(),
                         19)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(),
                         51)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-0.6841, -5.3454, 1.8254])
        self.assertEqual(prep.get_ligands()[18].get_identifier(), "15:3")

        # test write-out
        out_path = attach_root_path("tests/junk/corina_ligands_stereo.sdf")
        prep.write_ligands(path=out_path, format=_LP.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 93898)
