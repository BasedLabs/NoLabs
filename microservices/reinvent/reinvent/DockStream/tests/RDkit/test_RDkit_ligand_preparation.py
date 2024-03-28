import unittest
import os

from rdkit import Chem

from dockstream.core.ligand.ligand_input_parser import LigandInputParser
from dockstream.core.RDkit.RDkit_ligand_preparator import RDkitLigandPreparator
from dockstream.utils.enums.RDkit_enums import RDkitLigandPreparationEnum
from dockstream.utils.enums.rDock_enums import rDockDockingConfigurationEnum

from tests.tests_paths import PATHS_1UYD
from dockstream.utils.files_paths import attach_root_path
from dockstream.utils.smiles import read_smiles_file

_LP = RDkitLigandPreparationEnum()
_CE = rDockDockingConfigurationEnum()


class Test_RDkit_ligand_preparation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        self.ligands = []
        for mol in Chem.SDMolSupplier(attach_root_path(PATHS_1UYD.LIGANDS_SDF)):
            if mol is None:
                continue
            self.ligands.append(mol)
        self.ligands_smiles = list(read_smiles_file(attach_root_path(PATHS_1UYD.LIGANDS_SMILES_TXT), standardize=False))

    @classmethod
    def tearDownClass(cls):
        pass

    def test_coordinate_generation(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.TYPE: _LP.TYPE_RDKIT,
                _LP.PARAMS: {
                    _LP.EP_PARAMS_COORDGEN: {
                        _LP.EP_PARAMS_COORDGEN_METHOD: _LP.EP_PARAMS_COORDGEN_UFF,
                        _LP.EP_PARAMS_COORDGEN_UFF_MAXITERS: 350}
                }}
        lig_parser = LigandInputParser(smiles=self.ligands_smiles, **conf)
        prep = RDkitLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        # check 3D coordinate generation and subsequent alignment
        prep.generate3Dcoordinates(converged_only=False)
        self.assertEqual(prep.get_number_ligands(),
                         15)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(),
                         51)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-3.3083596377039033, -0.15090557789031064, 0.3769119664108664])

        # test write-out
        out_path = attach_root_path(("tests/junk/RDkit_ligands.sdf"))
        prep.write_ligands(path=out_path, format=_LP.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 62411)
        try:
            os.remove(out_path)
        except OSError:
            pass

    def test_aligning(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {
                    _LP.EP_PARAMS_COORDGEN: {
                        _LP.EP_PARAMS_COORDGEN_METHOD: _LP.EP_PARAMS_COORDGEN_UFF,
                        _LP.EP_PARAMS_COORDGEN_UFF_MAXITERS: 180}},
                    _LP.ALIGN: {_LP.ALIGN_MODE: _LP.ALIGN_MODE_INTERNAL,
                        _LP.ALIGN_REFERENCE_PATHS: [attach_root_path(PATHS_1UYD.LIGAND_PU8_SDF)],
                        _LP.ALIGN_REFERENCE_FORMAT: _LP.ALIGN_REFERENCE_FORMAT_SDF,
                        _LP.ALIGN_MINIMUM_SUBSTRUCTURE_RATIO: 0.2,
                        _LP.ALIGN_FAIL_ACTION: _LP.ALIGN_FAIL_DISCARD,
                        _LP.ALIGN_COMPLETE_RINGS_ONLY: True,
                        _LP.ALIGN_TETHERING: False}}
        lig_parser = LigandInputParser(smiles=self.ligands_smiles, **conf)
        prep = RDkitLigandPreparator(ligands=lig_parser.get_ligands(), **conf)
        prep.generate3Dcoordinates()

        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-3.3083596377039033, -0.15090557789031064, 0.3769119664108664])
        prep.align_ligands()
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [1.7116408325466044, 11.46692388076509, 27.921102089610418])

        # check other options (with and without Hs and convergence_only parameters)
        prep.generate3Dcoordinates(converged_only=False)
        self.assertEqual(len([True for lig in prep.get_ligands() if lig.get_molecule() is not None]),
                         15)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(),
                         51)
        self.assertEqual(prep.get_ligands()[2].get_molecule().GetNumAtoms(),
                         48)

        prep.generate3Dcoordinates(converged_only=True)
        self.assertEqual(len([True for lig in prep.get_ligands() if lig.get_molecule() is not None]),
                         11)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(),
                         51)
        self.assertTrue(prep.get_ligands()[2].get_molecule() is None)

    def test_tether_ligands(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {
                    _LP.EP_PARAMS_COORDGEN: {
                        _LP.EP_PARAMS_COORDGEN_METHOD: _LP.EP_PARAMS_COORDGEN_UFF,
                        _LP.EP_PARAMS_COORDGEN_UFF_MAXITERS: 350}},
                _LP.ALIGN: {_LP.ALIGN_MODE: _LP.ALIGN_MODE_INTERNAL,
                            _LP.ALIGN_REFERENCE_PATHS: [attach_root_path(PATHS_1UYD.LIGAND_PU8_PDB)],
                            _LP.ALIGN_REFERENCE_FORMAT: _LP.ALIGN_REFERENCE_FORMAT_PDB,
                            _LP.ALIGN_MINIMUM_SUBSTRUCTURE_RATIO: 0.2,
                            _LP.ALIGN_FAIL_ACTION: _LP.ALIGN_FAIL_DISCARD,
                            _LP.ALIGN_COMPLETE_RINGS_ONLY: True,
                            _LP.ALIGN_TETHERING: True}}
        lig_parser = LigandInputParser(smiles=self.ligands_smiles, **conf)
        prep = RDkitLigandPreparator(ligands=lig_parser.get_ligands(), **conf)
        prep.generate3Dcoordinates()
        prep.align_ligands()
        self.assertEqual(len(self.ligands), prep.get_number_ligands())
        self.assertEqual("1,2,3,4,5,6,7,8,9,10,11,12,13,16,19,17,18,14,15,20,21,27,26,25,24,22,23",
                         prep.get_ligands()[1].get_molecule().GetProp(_LP.TAG_ALIGNED_ATOMS))
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [1.7116408325466044, 11.46692388076509, 27.921102089610418])
