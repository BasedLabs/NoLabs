import tempfile
import unittest
import os
import rdkit.Chem as Chem

from dockstream.containers.docking_container import DockingContainer
from dockstream.core.Schrodinger.Glide_docker import Glide, AdvancedGlideKeywords, parse_maestro

from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum
from dockstream.utils.enums.docking_enum import ResultKeywordsEnum
from dockstream.utils.enums.Schrodinger_enums import SchrodingerDockingConfigurationEnum, SchrodingerOutputEnum, SchrodingerExecutablesEnum

from tests.tests_paths import PATHS_1UYD, PATH_SCHRODINGER_EXAMPLES
from dockstream.utils.files_paths import attach_root_path, any_in_file

from dockstream.core.ligand.ligand import Ligand
from dockstream.utils.smiles import to_smiles

from dockstream.utils.files_paths import lines_in_file

_LP = LigandPreparationEnum()
_CE = SchrodingerDockingConfigurationEnum()
_EE = SchrodingerExecutablesEnum()
_ROE = SchrodingerOutputEnum()
_RK = ResultKeywordsEnum()


# TODO: move "prefix_execution" to "config/tests_config/config.json"
class Test_Schrodinger_backend(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # specify absolute paths to the various input files
        cls.ligands_sdf_path = attach_root_path(PATHS_1UYD.LIGANDS_SDF_CORINA)
        cls.ligands_sdf_enumerated_path = attach_root_path(PATHS_1UYD.LIGANDS_WITH_ENUMERATION_SDF)
        cls.grid_file_path = attach_root_path(PATH_SCHRODINGER_EXAMPLES.GRIDFILE)

    def setUp(self):
        # "normal", embedded ligands
        list_ligands = []
        for lig_number, mol in enumerate(Chem.SDMolSupplier(self.ligands_sdf_path, removeHs=False)):
            list_ligands.append(Ligand(smile=to_smiles(mol),
                                       original_smile=to_smiles(mol),
                                       ligand_number=lig_number,
                                       enumeration=0,
                                       molecule=mol,
                                       mol_type=_LP.TYPE_RDKIT))
        self.ligands = list_ligands

        # enumerated, embedded ligands
        list_ligands = []
        for mol in Chem.SDMolSupplier(self.ligands_sdf_enumerated_path, removeHs=False):
            name = mol.GetProp("_Name").split(':')
            lig_number = int(name[0])
            enumeration = int(name[1])
            list_ligands.append(Ligand(smile=to_smiles(mol),
                                       original_smile=to_smiles(mol),
                                       ligand_number=lig_number,
                                       enumeration=enumeration,
                                       molecule=mol,
                                       mol_type=_LP.TYPE_RDKIT))
        self.ligands_enumerated = list_ligands

        # configuration
        self.conf = DockingContainer(conf={_CE.DOCKING: {_CE.DOCKING_RUNS: [
            {_CE.BACKEND: _CE.BACKEND_GLIDE,
             _CE.INPUT_POOLS: ["RDkit_pool1"],
             _CE.PARAMS: {
                 _CE.PARAMS_PREFIX_EXECUTION: "module load schrodinger/2019-4",
                 _CE.GLIDE_FLAGS: {
                     _EE.GLIDE_OVERWRITE: "",
                     _EE.GLIDE_NJOBS: 1,
                     _EE.GLIDE_HOST: "localhost",
                     _EE.GLIDE_TMPLAUNCHDIR: ""
                 },
                 _CE.GLIDE_KEYWORDS: {
                     _EE.GLIDE_AMIDE_MODE: "trans",
                     _EE.GLIDE_EXPANDED_SAMPLING: "True",
                     _EE.GLIDE_GRIDFILE: [self.grid_file_path],
                     _EE.GLIDE_NENHANCED_SAMPLING: '1',
                     _EE.GLIDE_POSE_OUTTYPE: _EE.GLIDE_POSE_OUTTYPE_LIGANDLIB,
                     _EE.GLIDE_POSES_PER_LIG: '3',
                     _EE.GLIDE_POSTDOCK_NPOSE: "25",
                     _EE.GLIDE_POSTDOCKSTRAIN: "True",
                     _EE.GLIDE_PRECISION: "SP",
                     _EE.GLIDE_REWARD_INTRA_HBONDS: "True"}}}]}})

    def test_Glide_input_file_preparation_simple(self):
        in_path = attach_root_path("tests/junk/glide_docking_simple.in")
        self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][0][_CE.PARAMS][_CE.GLIDE_KEYWORDS]["[CONSTRAINT_GROUP:1]"] = {
                    _EE.GLIDE_USE_CONS: "sulfonamide:1",
                    _EE.GLIDE_NREQUIRED_CONS: "ALL"}
        if os.path.isfile(in_path):
            os.remove(in_path)

        docking_run = 0
        run_parameters = self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][docking_run]
        docker = Glide(**run_parameters)
        keywords = docker._run_parameters[_CE.PARAMS][_CE.GLIDE_KEYWORDS]

        # test write-out with path specified
        docker._write_keywords_to_file(keywords=keywords, path=in_path)
        self.assertTrue(os.path.isfile(in_path))
        self.assertEqual(lines_in_file(in_path),
                         14)

        # test write-out without path specified (temporary file will be generated)
        tmp_in_path = docker._write_keywords_to_file(keywords=keywords)
        self.assertTrue(os.path.isfile(tmp_in_path))
        self.assertEqual(lines_in_file(tmp_in_path),
                         14)

    def test_Glide_input_file_preparation_constraint(self):
        in_path = attach_root_path("tests/junk/glide_docking_constraint.in")
        self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][0][_CE.PARAMS][_CE.GLIDE_KEYWORDS]["[CONSTRAINT_GROUP:1]"] = {
                    _EE.GLIDE_USE_CONS: "sulfonamide:1",
                    _EE.GLIDE_NREQUIRED_CONS: "ALL"}
        self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][0][_CE.PARAMS][_CE.GLIDE_KEYWORDS]["[FEATURE:1]"] = {
            "PATTERN1": "[N]#C 1 include",
            "PATTERN2": "[n] 1 include",
            "PATTERN3": "N(=N=N) 1 include",
            "PATTERN4": "N(=N)=N 1 include"
        }
        if os.path.isfile(in_path):
            os.remove(in_path)

        docking_run = 0
        run_parameters = self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][docking_run]
        docker = Glide(**run_parameters)
        keywords = docker._run_parameters[_CE.PARAMS][_CE.GLIDE_KEYWORDS]

        # test write-out with path specified
        docker._write_keywords_to_file(keywords=keywords, path=in_path)
        self.assertTrue(os.path.isfile(in_path))
        self.assertEqual(lines_in_file(in_path),
                         20)

        # test write-out without path specified (temporary file will be generated)
        tmp_in_path = docker._write_keywords_to_file(keywords=keywords)
        self.assertTrue(os.path.isfile(tmp_in_path))
        self.assertEqual(lines_in_file(tmp_in_path),
                         20)

    def test_Glide_docking_and_result_parsing(self):
        docking_run = 0
        run_parameters = self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][docking_run]
        docker = Glide(**run_parameters)
        docker.add_molecules(molecules=self.ligands)
        docker.dock()

        # test the docked ligand poses
        list_docked_ligands = docker.get_docked_ligands()
        del list_docked_ligands[0]
        list_docked_ligands = docker.get_docked_ligands()
        self.assertEqual(15, len(list_docked_ligands))

        self.assertListEqual([conf.GetProp(_ROE.GLIDE_DOCKING_SCORE) for conf in docker.get_docked_ligands()[1].get_conformers()],
                             ["-6.05413", "-5.68204", "-5.63135"])

        # placed here, as the input is the result dictionary
        df = docker.get_result()
        self.assertEqual(15, sum(df[_RK.DF_LOWEST_CONFORMER]))
        self.assertListEqual(list(df[_RK.DF_SCORE]),
                             [-8.29018, -7.37956, -7.18861, -6.05413, -5.68204, -5.63135, -7.58044, -7.40553, -7.08389,
                              -7.69184, -7.68062, -7.38606, -8.57887, -7.96752, -6.65871, -7.95027, -7.74089, -6.9786,
                              -8.31461, -7.7674, -7.19897, -7.68268, -7.66724, -7.11361, -8.36201, -7.56475, -7.32564,
                              -7.48451, -7.42463, -7.2352, -8.56583, -7.72547, -6.94599, -7.70743, -7.30964, -6.27659,
                              -8.59383, -8.16292, -7.07524, -8.79777, -8.66899, -6.72157, -8.36043, -8.11039, -7.79249])
        self.assertListEqual(list(df[_RK.DF_LOWEST_CONFORMER]),
                             [True, False, False, True, False, False, True, False, False, True, False, False, True,
                              False, False, True, False, False, True, False, False, True, False, False, True, False,
                              False, True, False, False, True, False, False, True, False, False, True, False, False,
                              True, False, False, True, False, False])

        # check the score extraction
        self.assertListEqual(docker.get_scores(best_only=True),
                             [-8.29018, -6.05413, -7.58044, -7.69184, -8.57887, -7.95027, -8.31461, -7.68268, -8.36201,
                              -7.48451, -8.56583, -7.70743, -8.59383, -8.79777, -8.36043])
        self.assertListEqual(docker.get_scores(best_only=False),
                             [-8.29018, -7.37956, -7.18861, -6.05413, -5.68204, -5.63135, -7.58044, -7.40553, -7.08389,
                              -7.69184, -7.68062, -7.38606, -8.57887, -7.96752, -6.65871, -7.95027, -7.74089, -6.9786,
                              -8.31461, -7.7674, -7.19897, -7.68268, -7.66724, -7.11361, -8.36201, -7.56475, -7.32564,
                              -7.48451, -7.42463, -7.2352, -8.56583, -7.72547, -6.94599, -7.70743, -7.30964, -6.27659,
                              -8.59383, -8.16292, -7.07524, -8.79777, -8.66899, -6.72157, -8.36043, -8.11039, -7.79249])

        # write out poses and check length
        out_path = attach_root_path("tests/junk/Glide_backend_docked_all.sdf")
        docker.write_docked_ligands(path=out_path, mode=_CE.OUTPUT_MODE_ALL)
        self.assertEqual(lines_in_file(out_path), 8514)
        out_path = attach_root_path("tests/junk/Glide_backend_docked_best_per_enumeration.sdf")
        docker.write_docked_ligands(path=out_path, mode=_CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertEqual(lines_in_file(out_path), 2838)

        # write out dataframe and check length
        out_path = attach_root_path("tests/junk/Glide_backend_docked_all.csv")
        docker.write_result(path=out_path, mode=_CE.OUTPUT_MODE_ALL)
        self.assertEqual(lines_in_file(out_path), 46)
        out_path = attach_root_path("tests/junk/Glide_backend_docked_best_per_enumeration.csv")
        docker.write_result(path=out_path, mode=_CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertEqual(lines_in_file(out_path), 16)

    def test_Glide_docking_and_result_parsing_parallelized(self):
        self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][0][_CE.PARAMS][_CE.PARALLELIZATION] = {
            _CE.PARALLELIZATION_NUMBER_CORES: 4}
        docking_run = 0
        run_parameters = self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][docking_run]
        docker = Glide(**run_parameters)
        docker.add_molecules(molecules=self.ligands)
        docker.dock()

        # test the docked ligand poses (make sure, the read-out is not affected by deletion, changes or additions)
        list_docked_ligands = docker.get_docked_ligands()
        del list_docked_ligands[0]
        list_docked_ligands = docker.get_docked_ligands()
        self.assertEqual(15, len(list_docked_ligands))

        self.assertListEqual([conf.GetProp(_ROE.GLIDE_DOCKING_SCORE) for conf in docker.get_docked_ligands()[1].get_conformers()],
                             ["-6.05413", "-5.68204", "-5.63135"])

        # placed here, as the input is the result dictionary
        df = docker.get_result()
        self.assertEqual(15, sum(df[_RK.DF_LOWEST_CONFORMER]))
        self.assertListEqual(list(df[_RK.DF_SCORE]),
                             [-8.29018, -7.37956, -7.18861, -6.05413, -5.68204, -5.63135, -7.58044, -7.40553, -7.08389,
                              -7.69184, -7.68062, -7.38606, -8.57887, -7.96752, -6.65871, -7.95027, -7.74089, -6.9786,
                              -8.31461, -7.7674, -7.19897, -7.68268, -7.66724, -7.11361, -8.36201, -7.56475, -7.32564,
                              -7.48451, -7.42463, -7.2352, -8.56583, -7.72547, -6.94599, -7.70743, -7.30964, -6.27659,
                              -8.59383, -8.16292, -7.07524, -8.79777, -8.66899, -6.72157, -8.36043, -8.11039, -7.79249])
        self.assertListEqual(list(df[_RK.DF_LOWEST_CONFORMER]),
                             [True, False, False, True, False, False, True, False, False, True, False, False, True,
                              False, False, True, False, False, True, False, False, True, False, False, True, False,
                              False, True, False, False, True, False, False, True, False, False, True, False, False,
                              True, False, False, True, False, False])

        # write out the result and check length
        out_path = attach_root_path("tests/junk/Glide_backend_docked_parallelized.sdf")
        docker.write_docked_ligands(path=out_path)
        self.assertEqual(lines_in_file(out_path), 8514)

        # write out dataframe and check length
        out_path = attach_root_path("tests/junk/Glide_backend_docked_parallelized.csv")
        docker.write_result(path=out_path)
        self.assertEqual(lines_in_file(out_path), 46)

    def test_Glide_docking_and_result_parsing_parallelized_subjob_restriction(self):
        self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][0][_CE.PARAMS][_CE.PARALLELIZATION] = {
            _CE.PARALLELIZATION_NUMBER_CORES: 4,
            _CE.PARALLELIZATION_MAXCOMPOUNDSPERSUBJOB: 1}
        docking_run = 0
        run_parameters = self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][docking_run]
        docker = Glide(**run_parameters)
        docker.add_molecules(molecules=self.ligands)
        docker.dock()

        # test the docked ligand poses (make sure, the read-out is not affected by deletion, changes or additions)
        list_docked_ligands = docker.get_docked_ligands()
        del list_docked_ligands[0]
        list_docked_ligands = docker.get_docked_ligands()
        self.assertEqual(15, len(list_docked_ligands))

        self.assertListEqual([conf.GetProp(_ROE.GLIDE_DOCKING_SCORE) for conf in docker.get_docked_ligands()[1].get_conformers()],
                             ["-6.05413", "-5.68204", "-5.63135"])

        # placed here, as the input is the result dictionary
        df = docker.get_result()
        self.assertEqual(15, sum(df[_RK.DF_LOWEST_CONFORMER]))
        self.assertListEqual(list(df[_RK.DF_SCORE]),
                             [-8.29018, -7.37956, -7.18861, -6.05413, -5.68204, -5.63135, -7.58044, -7.40553, -7.08389,
                              -7.69184, -7.68062, -7.38606, -8.57887, -7.96752, -6.65871, -7.95027, -7.74089, -6.9786,
                              -8.31461, -7.7674, -7.19897, -7.68268, -7.66724, -7.11361, -8.36201, -7.56475, -7.32564,
                              -7.48451, -7.42463, -7.2352, -8.56583, -7.72547, -6.94599, -7.70743, -7.30964, -6.27659,
                              -8.59383, -8.16292, -7.07524, -8.79777, -8.66899, -6.72157, -8.36043, -8.11039, -7.79249])
        self.assertListEqual(list(df[_RK.DF_LOWEST_CONFORMER]),
                             [True, False, False, True, False, False, True, False, False, True, False, False, True,
                              False, False, True, False, False, True, False, False, True, False, False, True, False,
                              False, True, False, False, True, False, False, True, False, False, True, False, False,
                              True, False, False, True, False, False])

        # write out the result and check length
        out_path = attach_root_path("tests/junk/Glide_backend_docked_parallelized.sdf")
        docker.write_docked_ligands(path=out_path)
        self.assertEqual(lines_in_file(out_path), 8514)

        # write out dataframe and check length
        out_path = attach_root_path("tests/junk/Glide_backend_docked_parallelized.csv")
        docker.write_result(path=out_path)
        self.assertEqual(lines_in_file(out_path), 46)

    def test_Glide_docking_and_result_parsing_parallelized_enumerated_ligands(self):
        self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][0][_CE.PARAMS][_CE.PARALLELIZATION] = {
            _CE.PARALLELIZATION_NUMBER_CORES: 4}
        docking_run = 0
        run_parameters = self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][docking_run]
        docker = Glide(**run_parameters)
        docker.add_molecules(molecules=self.ligands_enumerated)
        docker.dock()

        # test the docked ligand poses (make sure, the read-out is not affected by deletion, changes or additions)
        list_docked_ligands = docker.get_docked_ligands()
        del list_docked_ligands[0]
        list_docked_ligands = docker.get_docked_ligands()
        self.assertEqual(21, len(list_docked_ligands))

        self.assertListEqual([conf.GetProp(_ROE.GLIDE_DOCKING_SCORE) for conf in docker.get_docked_ligands()[16].get_conformers()],
                             ["-7.72026", "-7.41442", "-7.07512"])

        # placed here, as the input is the result dictionary; note that docking of the last 4 enumerations (one ligand)
        # failed, so only 17 (the 15 standard molecules and two enumerations of the additional one) are there; also, the
        # conformer lists of the last 4 enumerations are empty
        df = docker.get_result()
        self.assertEqual(17, sum(df[_RK.DF_LOWEST_CONFORMER]))
        self.assertListEqual(list(df[_RK.DF_SCORE]),
                             [-8.2844, -7.36993, -7.18515, -6.03692, -5.68207, -5.6232, -7.58758, -7.39131, -7.09023,
                              -7.68363, -7.66198, -6.80324, -8.4732, -7.68161, -6.38494, -7.95194, -7.74366, -6.98287,
                              -8.31415, -7.83831, -7.20183, -7.68931, -7.66631, -7.12742, -8.36191, -7.56978, -7.3317,
                              -7.48619, -7.43026, -7.23536, -9.47902, -9.06411, -8.67728, -9.74393, -9.0783, -8.86463,
                              -9.84685, -8.82573, -8.06811, -7.93591, -7.90722, -7.63321, -8.3517, -8.1201, -7.79452,
                              -6.95939, -6.85821, -6.8323, -7.72026, -7.41442, -7.07512])

        self.assertListEqual(list(df[_RK.DF_LOWEST_CONFORMER]),
                             [True, False, False, True, False, False, True, False, False, True, False, False, True,
                              False, False, True, False, False, True, False, False, True, False, False, True, False,
                              False, True, False, False, True, False, False, True, False, False, True, False, False,
                              True, False, False, True, False, False, True, False, False, True, False, False])

        # write out the result and check length
        out_path = attach_root_path("tests/junk/Glide_backend_docked_parallelized_enumerated_all.sdf")
        docker.write_docked_ligands(path=out_path, mode=_CE.OUTPUT_MODE_ALL)
        self.assertEqual(lines_in_file(out_path), 9462)
        out_path = attach_root_path("tests/junk/Glide_backend_docked_parallelized_enumerated_bestperenumeration.sdf")
        docker.write_docked_ligands(path=out_path, mode=_CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertEqual(lines_in_file(out_path), 3154)
        out_path = attach_root_path("tests/junk/Glide_backend_docked_parallelized_enumerated_bestperligand.sdf")
        docker.write_docked_ligands(path=out_path, mode=_CE.OUTPUT_MODE_BESTPERLIGAND)
        self.assertEqual(lines_in_file(out_path), 2999)

        # write out dataframes and check lengths
        out_path = attach_root_path("tests/junk/Glide_backend_docked_parallelized_enumerated_all.csv")
        docker.write_result(path=out_path, mode=_CE.OUTPUT_MODE_ALL)
        self.assertEqual(lines_in_file(out_path), 52)
        out_path = attach_root_path("tests/junk/Glide_backend_docked_parallelized_enumerated_bestperenumeration.csv")
        docker.write_result(path=out_path, mode=_CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertEqual(lines_in_file(out_path), 18)
        out_path = attach_root_path("tests/junk/Glide_backend_docked_parallelized_enumerated_bestperligand.csv")
        docker.write_result(path=out_path, mode=_CE.OUTPUT_MODE_BESTPERLIGAND)
        self.assertEqual(lines_in_file(out_path), 17)

    def test_Glide_docking_and_result_parsing_parallelized_constraints(self):
        self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][0][_CE.PARAMS][_CE.PARALLELIZATION] = {
            _CE.PARALLELIZATION_NUMBER_CORES: 4}
        self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][0][_CE.PARAMS][_CE.GLIDE_KEYWORDS][
            "[CONSTRAINT_GROUP:1]"] = {
            # use one constraint selectively, "A:SER:52:O(hbond)" would also be defined in this receptor grid
            _EE.GLIDE_USE_CONS: "A:ASP:93:OD1(hbond):1,",
            _EE.GLIDE_NREQUIRED_CONS: "ALL"}
        self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][0][_CE.PARAMS][_CE.GLIDE_KEYWORDS][
            "[FEATURE:1]"] = {
            "PATTERN1": "[#1][#7] 1 include"
        }

        self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][0][_CE.PARAMS][_CE.GLIDE_KEYWORDS][_EE.GLIDE_GRIDFILE] = [attach_root_path(PATH_SCHRODINGER_EXAMPLES.CONSTRAINTS_GRIDFILE)]
        docking_run = 0
        run_parameters = self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][docking_run]
        docker = Glide(**run_parameters)
        docker.add_molecules(molecules=self.ligands_enumerated)
        docker.dock()

        # test the docked ligand poses (make sure, the read-out is not affected by deletion, changes or additions)
        list_docked_ligands = docker.get_docked_ligands()
        del list_docked_ligands[0]
        list_docked_ligands = docker.get_docked_ligands()
        self.assertEqual(21, len(list_docked_ligands))

        # scores are slightly changed, because of
        # would be ["-7.41833", "-7.33341", "-7.32742"] without any activated constraints
        self.assertListEqual([conf.GetProp(_ROE.GLIDE_DOCKING_SCORE) for conf in docker.get_docked_ligands()[16].get_conformers()],
                             ["-7.57665", "-7.43678", "-6.54468"])

        # placed here, as the input is the result dictionary; note that docking of the last 4 enumerations (one ligand)
        # failed, so only 17 (the 15 standard molecules and two enumerations of the additional one) are there; also, the
        # conformer lists of the last 4 enumerations are empty
        df = docker.get_result()
        self.assertEqual(15, sum(df[_RK.DF_LOWEST_CONFORMER]))
        self.assertListEqual(list(df[_RK.DF_SCORE]),
                             [-7.52619, -7.45916, -5.28989, -6.50273, -4.79087, -7.25591, -7.24666, -6.94393, -7.53255,
                              -7.10787, -6.62053, -7.43503, -7.43478, -7.41845, -4.93922, -4.5561, -4.29516, -8.38331,
                              -8.17157, -8.06074, -7.43751, -7.19002, -7.18142, -8.0247, -7.47108, -8.00048, -7.70642,
                              -7.63619, -8.3103, -7.51347, -7.46712, -7.48016, -7.04096, -6.91424, -7.74463, -7.28038,
                              -7.09172, -6.49849, -6.39399, -6.19034, -7.57665, -7.43678, -6.54468])
        self.assertListEqual(list(df[_RK.DF_LOWEST_CONFORMER]),
                             [True, False, False, True, False, True, False, False, True, False, False, True, False,
                              False, True, False, False, True, False, False, True, False, False, True, False, True,
                              False, False, True, False, False, True, False, False, True, False, False, True, False,
                              False, True, False, False])

        # write out the result and check length
        out_path = attach_root_path("tests/junk/Glide_backend_docked_parallelized_enumerated_constraint_bestperligand.sdf")
        docker.write_docked_ligands(path=out_path, mode=_CE.OUTPUT_MODE_BESTPERLIGAND)
        self.assertEqual(lines_in_file(out_path), 2670)

        # write out dataframes and check lengths
        out_path = attach_root_path("tests/junk/Glide_backend_docked_parallelized_enumerated_constraint_all.csv")
        docker.write_result(path=out_path, mode=_CE.OUTPUT_MODE_ALL)
        self.assertEqual(lines_in_file(out_path), 44)

    def test_Glide_docking_with_token_guard(self):
        self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][0][_CE.PARAMS][_CE.GLIDE_TG] = {
            _CE.GLIDE_TG_PREFIX_EXECUTION: "module load schrodinger/2019-4",
            _CE.GLIDE_TG_TOKEN_POOLS: {
                "GLIDE_SP_DOCKING": 4
            },
            _CE.GLIDE_TG_WAIT_INTERVAL: 10,
            _CE.GLIDE_TG_WAIT_LIMIT: 60}
        docking_run = 0
        run_parameters = self.conf[_CE.DOCKING][_CE.DOCKING_RUNS][docking_run]
        docker = Glide(**run_parameters)
        docker.add_molecules(molecules=self.ligands)
        docker.dock()

        # test the docked ligand poses (make sure, the read-out is not affected by deletion, changes or additions)
        list_docked_ligands = docker.get_docked_ligands()
        del list_docked_ligands[0]
        list_docked_ligands = docker.get_docked_ligands()
        self.assertEqual(15, len(list_docked_ligands))

        self.assertListEqual([conf.GetProp(_ROE.GLIDE_DOCKING_SCORE) for conf in docker.get_docked_ligands()[1].get_conformers()],
                             ["-6.05413", "-5.68204", "-5.63135"])

        # placed here, as the input is the result dictionary
        df = docker.get_result()
        self.assertEqual(15, sum(df[_RK.DF_LOWEST_CONFORMER]))
        self.assertListEqual(list(df[_RK.DF_SCORE]),
                             [-8.29018, -7.37956, -7.18861, -6.05413, -5.68204, -5.63135, -7.58044, -7.40553, -7.08389,
                              -7.69184, -7.68062, -7.38606, -8.57887, -7.96752, -6.65871, -7.95027, -7.74089, -6.9786,
                              -8.31461, -7.7674, -7.19897, -7.68268, -7.66724, -7.11361, -8.36201, -7.56475, -7.32564,
                              -7.48451, -7.42463, -7.2352, -8.56583, -7.72547, -6.94599, -7.70743, -7.30964, -6.27659,
                              -8.59383, -8.16292, -7.07524, -8.79777, -8.66899, -6.72157, -8.36043, -8.11039, -7.79249])
        self.assertListEqual(list(df[_RK.DF_LOWEST_CONFORMER]),
                             [True, False, False, True, False, False, True, False, False, True, False, False, True,
                              False, False, True, False, False, True, False, False, True, False, False, True, False,
                              False, True, False, False, True, False, False, True, False, False, True, False, False,
                              True, False, False, True, False, False])

        # check the score extraction
        self.assertListEqual(docker.get_scores(best_only=True),
                             [-8.29018, -6.05413, -7.58044, -7.69184, -8.57887, -7.95027, -8.31461, -7.68268, -8.36201,
                              -7.48451, -8.56583, -7.70743, -8.59383, -8.79777, -8.36043])
        self.assertListEqual(docker.get_scores(best_only=False),
                             [-8.29018, -7.37956, -7.18861, -6.05413, -5.68204, -5.63135, -7.58044, -7.40553, -7.08389,
                              -7.69184, -7.68062, -7.38606, -8.57887, -7.96752, -6.65871, -7.95027, -7.74089, -6.9786,
                              -8.31461, -7.7674, -7.19897, -7.68268, -7.66724, -7.11361, -8.36201, -7.56475, -7.32564,
                              -7.48451, -7.42463, -7.2352, -8.56583, -7.72547, -6.94599, -7.70743, -7.30964, -6.27659,
                              -8.59383, -8.16292, -7.07524, -8.79777, -8.66899, -6.72157, -8.36043, -8.11039, -7.79249])

        # write out poses and check length
        out_path = attach_root_path("tests/junk/Glide_backend_docked.sdf")
        docker.write_docked_ligands(path=out_path)
        self.assertEqual(lines_in_file(out_path), 8514)

        # write out dataframe and check length
        out_path = attach_root_path("tests/junk/Glide_backend_docked.csv")
        docker.write_result(path=out_path)
        self.assertEqual(lines_in_file(out_path), 46)
