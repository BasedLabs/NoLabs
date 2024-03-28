import unittest
import os

from dockstream.core.OpenEye.OpenEye_ligand_preparator import OpenEyeLigandPreparator
from dockstream.utils.enums.OpenEye_enums import OpenEyeDockingConfigurationEnum, OpenEyeLigandPreparationEnum

from dockstream.core.ligand.ligand_input_parser import LigandInputParser

from tests.tests_paths import PATHS_1UYD, MAIN_CONFIG
from dockstream.utils.files_paths import attach_root_path
from dockstream.utils.smiles import read_smiles_file

_LP = OpenEyeLigandPreparationEnum()
_CE = OpenEyeDockingConfigurationEnum()


class Test_OpenEye_ligand_preparation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if "OE_LICENSE" in MAIN_CONFIG:
            os.environ["OE_LICENSE"] = MAIN_CONFIG["OE_LICENSE"]

    def setUp(self):
        self.ligands_smiles = list(read_smiles_file(attach_root_path(PATHS_1UYD.LIGANDS_SMILES_TXT),
                                                    standardize=False))

    @classmethod
    def tearDownClass(cls):
        pass

    def test_initializations(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.TYPE: _LP.TYPE_OPENEYE,
                _LP.PARAMS: {}}

        # initialize using a list of (RDkit) Mol objects
        lig_parser = LigandInputParser(smiles=self.ligands_smiles, **conf)
        prep = OpenEyeLigandPreparator(ligands=lig_parser.get_ligands(), **conf)
        self.assertEqual(15, prep.get_number_ligands())

        # initialize using a list of SMILES
        prep = OpenEyeLigandPreparator(ligands=lig_parser.get_ligands(), **conf)
        self.assertEqual(15, prep.get_number_ligands())

    def test_coordinate_generation(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {}}
        lig_parser = LigandInputParser(smiles=self.ligands_smiles, **conf)
        prep = OpenEyeLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        # check 3D coordinate generation
        prep.generate3Dcoordinates()
        self.assertEqual(15, prep.get_number_ligands())
        self.assertEqual((-3.754112720489502, -2.629582166671753, 6.873664855957031),
                         prep.get_ligands()[1].get_molecule().GetCoords()[0])

    def test_coordinate_generation_fails(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {}}
        test = self.ligands_smiles
        test.append("abc")
        lig_parser = LigandInputParser(smiles=test, **conf)
        prep = OpenEyeLigandPreparator(ligands=lig_parser.get_ligands(), **conf)
        self.assertEqual(16, prep.get_number_ligands())

        # check 3D coordinate generation
        prep.generate3Dcoordinates()
        self.assertEqual(16, prep.get_number_ligands())
        self.assertEqual((-3.754112720489502, -2.629582166671753, 6.873664855957031),
                         prep.get_ligands()[1].get_molecule().GetCoords()[0])
        self.assertEqual(len([True for lig in prep.get_ligands() if lig.get_molecule() is not None]),
                         15)

    def test_align_ligands(self):
        # note, that OpenEye is also able to align the molecules to a reference ligand, but this is stored in the
        # receptor, so no alignment needs to be done; here, the internal RDkit-based alignment is performed
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {},
                _LP.ALIGN: {_LP.ALIGN_MODE: _LP.ALIGN_MODE_INTERNAL,
                    _LP.ALIGN_REFERENCE_PATHS: [attach_root_path(PATHS_1UYD.LIGAND_PU8_SDF)],
                    _LP.ALIGN_REFERENCE_FORMAT: _LP.ALIGN_REFERENCE_FORMAT_SDF,
                    _LP.ALIGN_MINIMUM_SUBSTRUCTURE_RATIO: 0.2,
                    _LP.ALIGN_FAIL_ACTION: _LP.ALIGN_FAIL_DISCARD,
                    _LP.ALIGN_COMPLETE_RINGS_ONLY: True}}
        lig_parser = LigandInputParser(smiles=self.ligands_smiles, **conf)
        prep = OpenEyeLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        # generate coordinates
        prep.generate3Dcoordinates()

        # align ligands
        self.assertEqual((-5.948095798492432, 1.2310136556625366, -6.726933479309082),
                         prep.get_ligands()[0].get_molecule().GetCoords()[0])
        prep.align_ligands()
        self.assertEqual((0.8390467762947083, 6.771180152893066, 19.497631072998047),
                         prep.get_ligands()[0].get_molecule().GetCoords()[0])
