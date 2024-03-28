import unittest
import os

from dockstream.core.Schrodinger.Ligprep_ligand_preparator import LigprepLigandPreparator

from dockstream.core.ligand.ligand_input_parser import LigandInputParser

from dockstream.utils.enums.stereo_enumeration_enums import StereoEnumerationEnum
from dockstream.utils.enums.Schrodinger_enums import LigprepLigandPreparationEnum

from tests.tests_paths import PATHS_1UYD
from dockstream.utils.files_paths import attach_root_path
from dockstream.utils.smiles import read_smiles_file

_LP = LigprepLigandPreparationEnum()
_TE = StereoEnumerationEnum()


# TODO: move "prefix_execution" to "config/tests_config/config.json"
class Test_Ligprep_ligand_preparation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        self.smiles = list(read_smiles_file(attach_root_path(PATHS_1UYD.LIGANDS_SMILES_SMI_TAUT_ENUM),
                                            standardize=False))

    @classmethod
    def tearDownClass(cls):
        pass

    def test_coordinate_generation(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {
                    _LP.PREFIX_EXECUTION: "module load schrodinger/2020-4"}}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = LigprepLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        prep.generate3Dcoordinates()
        self.assertEqual(prep.get_number_ligands(),
                         19)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(),
                         51)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-7.3614, 1.1703, -3.0363])

        # check, that all could be embedded successfully
        self.assertEqual(len([True for lig in prep.get_ligands() if lig.get_molecule() is not None]),
                         18)

        # test write-out
        out_path = attach_root_path("tests/junk/ligands_ligprep.sdf")
        prep.write_ligands(path=out_path, format=_LP.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 83000)

    def test_coordinate_generation_AWS(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {
                    _LP.COMMAND_LINE_PARAMETERS: {"-HOST": "cpu-only"},
                    _LP.PREFIX_EXECUTION: "module load schrodinger/2020-4-js-aws"}}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = LigprepLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        prep.generate3Dcoordinates()
        self.assertEqual(prep.get_number_ligands(),
                         19)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(),
                         51)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-7.3614, 1.1703, -3.0363])

        # check, that all could be embedded successfully
        self.assertEqual(len([True for lig in prep.get_ligands() if lig.get_molecule() is not None]),
                         18)

        # test write-out
        out_path = attach_root_path("tests/junk/ligands_ligprep_aws.sdf")
        prep.write_ligands(path=out_path, format=_LP.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 83000)

    def test_failed_generation(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {
                    _LP.PREFIX_EXECUTION: "module load schrodinger/2019-4"}}
        lig_parser = LigandInputParser(smiles=['C#CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2Cl)nc2c(N)ncnc21', "abc", "another"],
                                       **conf)
        prep = LigprepLigandPreparator(ligands=lig_parser.get_ligands(), **conf)
        prep.generate3Dcoordinates()
        self.assertEqual(len([True for lig in prep.get_ligands() if lig.get_molecule() is not None]),
                         1)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(),
                         51)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-7.3614, 1.1703, -3.0363])

    def test_coordinate_generation_parallelized(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {
                    _LP.PARALLELIZATION: {
                        _LP.PARALLELIZATION_NUMBER_CORES: 4
                    },
                    _LP.PREFIX_EXECUTION: "module load schrodinger/2019-4"}}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = LigprepLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        prep.generate3Dcoordinates()
        self.assertEqual(prep.get_number_ligands(),
                         19)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(),
                         51)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-7.3614, 1.1703, -3.0363])

        # test write-out
        out_path = attach_root_path("tests/junk/ligands_ligprep_parallelized.sdf")
        prep.write_ligands(path=out_path, format=_LP.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 83000)

    def test_coordinate_generation_with_epik(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {
                    _LP.PREFIX_EXECUTION: "module load schrodinger/2019-4",
                    _LP.FORCE_FIELD: _LP.FORCE_FIELD_OPLS3e,
                    _LP.USE_EPIK: {
                        _LP.TARGET_PH: 7.0,
                        _LP.PH_TOLERANCE: 2.0
                    }}}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = LigprepLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        prep.generate3Dcoordinates()
        self.assertEqual(prep.get_number_ligands(),
                         20)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(),
                         51)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-6.2899, -1.2605, -0.4723])

        # test write-out
        out_path = attach_root_path("tests/junk/ligands_ligprep_epik.sdf")
        prep.write_ligands(path=out_path, format=_LP.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 100000)

    def test_coordinate_generation_with_epik_parameters_set(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {
                    _LP.PREFIX_EXECUTION: "module load schrodinger/2019-4",
                    _LP.FORCE_FIELD: _LP.FORCE_FIELD_OPLS3e,
                    _LP.USE_EPIK: {
                        _LP.TARGET_PH: 4.0,
                        _LP.PH_TOLERANCE: 1.0
                    }}}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = LigprepLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        prep.generate3Dcoordinates()
        self.assertEqual(prep.get_number_ligands(),
                         30)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(),
                         52)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [0.189, -1.0565, 4.9592])

    def test_coordinate_generation_with_epik_and_filtering(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {
                    _LP.PREFIX_EXECUTION: "module load schrodinger/2019-4",
                    _LP.FORCE_FIELD: _LP.FORCE_FIELD_OPLS3e,
                    _LP.USE_EPIK: {
                        _LP.TARGET_PH: 7.0,
                        _LP.PH_TOLERANCE: 2.0
                    },
                    _LP.FILTER_FILE: {
                        "Total_charge": "!= 0"}}}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = LigprepLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        prep.generate3Dcoordinates()
        self.assertListEqual([lig.get_smile() for lig in prep.get_ligands()],
                             ['[H]C#CC([H])([H])C([H])([H])C([H])([H])n1c(C([H])([H])c2c([H])c(OC([H])([H])[H])c(OC([H])([H])[H])c(OC([H])([H])[H])c2Cl)nc2c(N([H])[H])nc([H])nc21',
                              'CXXC',
                              '[H]c1nc(N([H])[H])c2nc(C([H])([H])c3c([H])c(OC([H])([H])[H])c(OC([H])([H])[H])c(OC([H])([H])[H])c3[H])n(C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])c2n1',
                              '[H]c1nc(N([H])[H])c2nc(C([H])([H])c3c([H])c(OC([H])([H])[H])c([H])c([H])c3OC([H])([H])[H])n(C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])c2n1',
                              '[H]c1nc(N([H])[H])c2nc(C([H])([H])c3c([H])c([H])c([H])c(OC([H])([H])[H])c3[H])n(C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])c2n1',
                              '[H]C#CC([H])([H])C([H])([H])C([H])([H])n1c(C([H])([H])c2c([H])c(OC([H])([H])[H])c(OC([H])([H])[H])c(OC([H])([H])[H])c2Cl)nc2c(N([H])[H])nc(F)nc21',
                              '[H]c1nc(N([H])[H])c2nc(C([H])([H])c3c([H])c([H])c(OC([H])([H])[H])c([H])c3[H])n(C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])c2n1',
                              '[H]c1nc(N([H])[H])c2nc(C([H])([H])c3c([H])c([H])c4c(c3[H])OC([H])([H])O4)n(C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])c2n1',
                              '[H]c1c([H])c(OC([H])([H])[H])c(C([H])([H])c2nc3c(N([H])[H])nc(F)nc3n2C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])c([H])c1OC([H])([H])[H]',
                              '[H]c1c([H])c(C([H])([H])c2nc3c(N([H])[H])nc(F)nc3n2C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])c([H])c2c1OC([H])([H])O2',
                              '[H]C#CC([H])([H])C([H])([H])C([H])([H])n1c(C([H])([H])c2c([H])c(OC([H])([H])[H])c([H])c([H])c2OC([H])([H])[H])nc2c(N([H])[H])nc(F)nc21',
                              'CC(C)NCCCn1c(Cc2cc3c(cc2I)OCO3)nc2c(N)nc(F)nc21', # these three have been filtered out, as they have +1 total charge
                              'CC(C)NCCCn1c(Sc2cc3c(cc2Br)OCO3)nc2c(N)ncnc21',
                              'CC(C)NCCCn1c(Sc2cc3c(cc2I)OCO3)nc2c(N)ncnc21',
                              '[H]c1c([H])c(OC([H])([H])[H])c(C([H])([H])c2nc3c(N([H])[H])nc(F)nc3n2[H])c([H])c1OC([H])([H])[H]',
                              '[H]c1c([H])c(OC([H])([H])[H])c(C([H])([H])c2nc3nc(F)nc(N([H])[H])c3n2[H])c([H])c1OC([H])([H])[H]',
                              '[H]c1nc(N([H])[H])c2nc(C([H])([H])c3c([H])c4c(c([H])c3Br)OC([H])([H])O4)c(N([H])C([H])([H])c3c([H])c([H])c([H])c([H])c3[H])n2c1[H]',
                              '[H]c1sc(N=C(N([H])[H])N([H])[H])nc1-c1c([H])c(C([H])([H])[H])n(C([H])([H])[H])c1[H]',
                              '[H]c1c([H])c([H])c2c(c1[H])C(=O)C([H])(C1=NS(=O)(=O)c3c([H])c(N([H])S(=O)(=O)C([H])([H])[H])c([H])c([H])c3N1[H])C(=O)C2(N([H])C([H])([H])c1c([H])c([H])c([H])c(OC([H])([H])[H])c1[H])C([H])([H])C([H])([H])C(C([H])([H])[H])(C([H])([H])[H])C([H])([H])[H]'])
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(),
                         51)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-6.2899, -1.2605, -0.4723])

        # test write-out
        out_path = attach_root_path("tests/junk/ligands_ligprep_epik_and_filtering.sdf")
        prep.write_ligands(path=out_path, format=_LP.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertGreater(100000, stat_inf.st_size)
        self.assertGreater(stat_inf.st_size, 60000)

    def test_coordinate_generation_with_Ligprep_stereo(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {
                    _LP.PREFIX_EXECUTION: "module load schrodinger/2019-4",
                    _LP.CHIRALITY: {}}}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = LigprepLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        prep.generate3Dcoordinates()
        self.assertEqual(prep.get_number_ligands(),
                         22)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(),
                         51)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-7.3614, 1.1703, -3.0363])

        # test write-out
        out_path = attach_root_path("tests/junk/ligands_ligprep_chiral.sdf")
        prep.write_ligands(path=out_path, format=_LP.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 106000)

    def test_coordinate_generation_with_Ligprep_stereo_command_line_parameters(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {
                    _LP.PREFIX_EXECUTION: "module load schrodinger/2019-4",
                    _LP.COMMAND_LINE_PARAMETERS: {"-g": ""},
                    _LP.CHIRALITY: {}}}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = LigprepLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        prep.generate3Dcoordinates()
        self.assertEqual(prep.get_number_ligands(),
                         19)
        self.assertEqual(prep.get_ligands()[18].get_molecule().GetNumAtoms(),
                         81)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-7.3614, 1.1703, -3.0363])

        # test write-out
        out_path = attach_root_path("tests/junk/ligands_ligprep_chiral_CLI_parameters_test.sdf")
        prep.write_ligands(path=out_path, format=_LP.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 83000)

    def test_coordinate_generation_with_RDkit_stereo(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST,
                            _TE.STEREO_ENUM: {
                                _TE.STEREO_ENUM_BACKEND: _TE.STEREO_ENUM_BACKEND_RDKIT}
                            },
                _LP.PARAMS: {
                    _LP.PREFIX_EXECUTION: "module load schrodinger/2019-4"}}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = LigprepLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        prep.generate3Dcoordinates()
        self.assertEqual(prep.get_number_ligands(),
                         34)
        self.assertEqual(prep.get_ligands()[0].get_molecule().GetNumAtoms(),
                         51)
        self.assertListEqual(list(prep.get_ligands()[33].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-4.5416, 1.1427, 2.7519])

        # test write-out
        out_path = attach_root_path("tests/junk/ligands_ligprep_rdkit_stereo.sdf")
        prep.write_ligands(path=out_path, format=_LP.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 170000)

    def test_aligning_keep(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {
                    _LP.PREFIX_EXECUTION: "module load schrodinger/2019-4",
                    _LP.USE_EPIK: {
                        _LP.TARGET_PH: 7.0,
                        _LP.PH_TOLERANCE: 2.0
                    }
                },
                _LP.ALIGN: {_LP.ALIGN_MODE: _LP.ALIGN_MODE_INTERNAL,
                                 _LP.ALIGN_REFERENCE_PATHS: [attach_root_path(PATHS_1UYD.LIGAND_PU8_SDF)],
                                 _LP.ALIGN_REFERENCE_FORMAT: _LP.ALIGN_REFERENCE_FORMAT_SDF,
                                 _LP.ALIGN_MINIMUM_SUBSTRUCTURE_RATIO: 0.2,
                                 _LP.ALIGN_FAIL_ACTION: _LP.ALIGN_FAIL_KEEP,
                                 _LP.ALIGN_COMPLETE_RINGS_ONLY: True}}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = LigprepLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        prep.generate3Dcoordinates()

        # entry number 15 (the 16th) will fail to align, but as ALIGN_FAIL_ACTION is KEEP, it should still have its
        # original coordinates
        # note: the z-values close to 0 for entry number 15 is fine (very small planar molecule)
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-7.3614, 1.1703, -3.0363])
        self.assertListEqual(list(prep.get_ligands()[17].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-3.5283, 1.3584, 0.0])
        self.assertEqual(len(prep.get_ligands()), 20)
        prep.align_ligands()
        self.assertEqual(len(prep.get_ligands()), 20)
        self.assertListEqual(list(prep.get_ligands()[17].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-3.5283, 1.3584, 0.0])
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-0.7020060673106414, 4.754538082396168, 21.483911269630926])

        # test write-out
        out_path = attach_root_path("tests/junk/ligands_ligprep_align_keep.sdf")
        prep.write_ligands(path=out_path, format=_LP.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 102000)

    def test_aligning_discard(self):
        conf = {_LP.POOLID: "testPool",
                _LP.INPUT: {_LP.INPUT_TYPE: _LP.INPUT_TYPE_LIST},
                _LP.PARAMS: {
                    _LP.PREFIX_EXECUTION: "module load schrodinger/2019-4",
                    _LP.USE_EPIK: {
                        _LP.TARGET_PH: 7.0,
                        _LP.PH_TOLERANCE: 2.0
                    }},
                _LP.ALIGN: {_LP.ALIGN_MODE: _LP.ALIGN_MODE_INTERNAL,
                                 _LP.ALIGN_REFERENCE_PATHS: [attach_root_path(PATHS_1UYD.LIGAND_PU8_SDF)],
                                 _LP.ALIGN_REFERENCE_FORMAT: _LP.ALIGN_REFERENCE_FORMAT_SDF,
                                 _LP.ALIGN_MINIMUM_SUBSTRUCTURE_RATIO: 0.2,
                                 _LP.ALIGN_FAIL_ACTION: _LP.ALIGN_FAIL_DISCARD,
                                 _LP.ALIGN_COMPLETE_RINGS_ONLY: True}}
        lig_parser = LigandInputParser(smiles=self.smiles, **conf)
        prep = LigprepLigandPreparator(ligands=lig_parser.get_ligands(), **conf)

        prep.generate3Dcoordinates()

        # entry number 15 (the 16th) will fail to align, as as ALIGN_FAIL_ACTION is DISCARD the "molecule" is "None"
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-7.3614, 1.1703, -3.0363])
        self.assertListEqual(list(prep.get_ligands()[17].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-3.5283, 1.3584, 0.0])
        self.assertEqual(len(prep.get_ligands()), 20)
        prep.align_ligands()
        self.assertEqual(len(prep.get_ligands()), 20)
        self.assertIsNone(prep.get_ligands()[17].get_molecule())
        self.assertListEqual(list(prep.get_ligands()[0].get_molecule().GetConformer(0).GetPositions()[0]),
                             [-0.7020060673106414, 4.754538082396168, 21.483911269630926])

        # test write-out
        out_path = attach_root_path("tests/junk/ligands_ligprep_align_discard.sdf")
        prep.write_ligands(path=out_path, format=_LP.OUTPUT_FORMAT_SDF)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 94000)
