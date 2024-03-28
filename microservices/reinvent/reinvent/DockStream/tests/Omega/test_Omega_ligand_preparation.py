import unittest
import os

from dockstream.core.OpenEyeHybrid.Omega_ligand_preparator import OmegaLigandPreparator

from dockstream.core.ligand.ligand_input_parser import LigandInputParser

from dockstream.utils.enums.docking_enum import DockingConfigurationEnum
from dockstream.utils.enums.Omega_enums import OmegaExecutablesEnum

from tests.tests_paths import PATHS_1UYD
from dockstream.utils.files_paths import attach_root_path
from dockstream.utils.smiles import read_smiles_file

_DE = DockingConfigurationEnum()
_OE = OmegaExecutablesEnum()


# TODO: move "prefix_execution" to "config/tests_config/config.json"
# TODO: add macrocycle mode unit test
class Test_Omega_ligand_preparation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        self.smiles = list(read_smiles_file(attach_root_path(PATHS_1UYD.LIGANDS_SMILES_TXT), standardize=False))

    @classmethod
    def tearDownClass(cls):
        pass

    def test_classic_coordinate_generation(self):
        conf = {_OE.POOLID: "testPool",
                _OE.INPUT: {_OE.INPUT_TYPE: _OE.INPUT_TYPE_LIST},
                _OE.PARAMS: {_OE.PREFIX_EXECUTION: "module load omega",
                             _OE.MODE: _OE.CLASSIC}}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = OmegaLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        # check 3D coordinate generation
        prep.generate3Dcoordinates()
        self.assertEqual(prep.get_number_ligands(), 15)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(), 51)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-5.9481, 1.231, -6.7269])

        # test write-out
        out_path = attach_root_path("tests/junk/omega_classic_ligands.sdf")
        prep.write_ligands(path=out_path, format=_OE.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 63824)

    def test_classic_coordinate_generation_fails(self):
        conf = {_OE.POOLID: "testPool",
                _OE.INPUT: {_OE.INPUT_TYPE: _OE.INPUT_TYPE_LIST},
                _OE.PARAMS: {_OE.PREFIX_EXECUTION: "module load omega",
                             _OE.MODE: _OE.CLASSIC}}
        test = self.smiles + ["abc"]
        lig_parser = LigandInputParser(smiles=test, **conf)
        prep = OmegaLigandPreparator(ligands=lig_parser.get_ligands(), **conf)
        self.assertEqual(prep.get_number_ligands(), 16)

        # check 3D coordinate generation
        prep.generate3Dcoordinates()
        self.assertEqual(prep.get_number_ligands(), 16)
        self.assertEqual(len([True for lig in prep.get_ligands() if lig.get_molecule() is not None]),
                         15)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(), 51)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-5.9481, 1.231, -6.7269])

        # test write-out
        out_path = attach_root_path("tests/junk/omega_classic_ligands.sdf")
        prep.write_ligands(path=out_path, format=_OE.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 63824)

    def test_classic_coordinate_generation_parallelized(self):
        conf = {_OE.POOLID: "testPool",
                _OE.INPUT: {_OE.INPUT_TYPE: _OE.INPUT_TYPE_LIST},
                _OE.PARAMS: {
                    _OE.PARALLELIZATION: {
                        _OE.PARALLELIZATION_NUMBER_CORES: 4
                    },
                    _OE.PREFIX_EXECUTION: "module load omega",
                    _OE.MODE: _OE.CLASSIC}}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = OmegaLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        # check 3D coordinate generation
        prep.generate3Dcoordinates()
        self.assertEqual(prep.get_number_ligands(), 15)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(), 51)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-5.9481, 1.231, -6.7269])

        # test write-out
        out_path = attach_root_path("tests/junk/omega_classic_parallelized_ligands.sdf")
        prep.write_ligands(path=out_path, format=_OE.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 63824)

    def test_rocs_coordinate_generation(self):
        conf = {_OE.POOLID: "testPool",
                _OE.INPUT: {_OE.INPUT_TYPE: _OE.INPUT_TYPE_LIST},
                _OE.PARAMS: {_OE.PREFIX_EXECUTION: "module load omega",
                             _OE.MODE: _OE.ROCS}}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = OmegaLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        # check 3D coordinate generation
        prep.generate3Dcoordinates()
        self.assertEqual(prep.get_number_ligands(), 15)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(), 51)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-5.9481, 1.231, -6.7269])

        # test write-out
        out_path = attach_root_path("tests/junk/omega_rocs_ligands.sdf")
        prep.write_ligands(path=out_path, format=_OE.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 63824)

    def test_pose_coordinate_generation(self):
        conf = {_OE.POOLID: "testPool",
                _OE.INPUT: {_OE.INPUT_TYPE: _OE.INPUT_TYPE_LIST},
                _OE.PARAMS: {_OE.PREFIX_EXECUTION: "module load omega",
                             _OE.MODE: _OE.POSE}}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = OmegaLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        # check 3D coordinate generation
        prep.generate3Dcoordinates()
        self.assertEqual(prep.get_number_ligands(), 15)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(), 51)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-5.9481, 1.231, -6.7269])

        # test write-out
        out_path = attach_root_path("tests/junk/omega_pose_ligands.sdf")
        prep.write_ligands(path=out_path, format=_OE.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 63824)

    def test_dense_coordinate_generation(self):
        conf = {_OE.POOLID: "testPool",
                _OE.INPUT: {_OE.INPUT_TYPE: _OE.INPUT_TYPE_LIST},
                _OE.PARAMS: {_OE.PREFIX_EXECUTION: "module load omega",
                             _OE.MODE: _OE.DENSE}}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = OmegaLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        # check 3D coordinate generation
        prep.generate3Dcoordinates()
        self.assertEqual(prep.get_number_ligands(), 15)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(), 51)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-5.9481, 1.231, -6.7269])

        # test write-out
        out_path = attach_root_path("tests/junk/omega_dense_ligands.sdf")
        prep.write_ligands(path=out_path, format=_OE.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 63824)
