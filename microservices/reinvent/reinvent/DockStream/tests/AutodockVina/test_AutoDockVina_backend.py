import unittest
import os
import rdkit.Chem as Chem

from dockstream.core.AutodockVina.AutodockVina_docker import AutodockVina, AutodockVinaParameters, SearchSpace
from dockstream.core.Schrodinger.Glide_docker import Parallelization

from dockstream.utils.enums.AutodockVina_enums import AutodockVinaDockingConfigurationEnum, \
                                                  AutodockResultKeywordsEnum
from dockstream.utils.enums.RDkit_enums import RDkitLigandPreparationEnum

from tests.tests_paths import PATHS_1UYD, PATH_AUTODOCKVINA_EXAMPLES
from dockstream.utils.files_paths import attach_root_path
from dockstream.core.ligand.ligand import Ligand
from dockstream.utils.smiles import to_smiles

from dockstream.utils.files_paths import lines_in_file


class Test_AutoDockVina_backend(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._CE = AutodockVinaDockingConfigurationEnum()
        cls._LP = RDkitLigandPreparationEnum()
        cls._ROE = AutodockResultKeywordsEnum()

        # specify absolute paths to the various input files
        cls.receptor_path = attach_root_path(PATH_AUTODOCKVINA_EXAMPLES.RECEPTOR)
        cls.ligands_with_hydrogens_sdf = attach_root_path(PATHS_1UYD.LIGANDS_SDF_WITH_HYDROGENS)
        cls.ligand_enum_sdf = attach_root_path(PATHS_1UYD.LIGANDS_WITH_ENUMERATION_SDF)

        cls._folder_dir = attach_root_path(PATH_AUTODOCKVINA_EXAMPLES.BACKEND_TESTS_FOLDER)
        if not os.path.isdir(cls._folder_dir):
            os.makedirs(cls._folder_dir)

    def setUp(self):
        list_ligands_with_hydrogens = []
        # Populate a list with the raw information in the input ligands.sdf file
        for lig_number, mol in enumerate(Chem.SDMolSupplier(self.ligands_with_hydrogens_sdf, removeHs=False)):
            list_ligands_with_hydrogens.append(Ligand(smile=to_smiles(mol),
                                                      original_smile=to_smiles(mol),
                                                      ligand_number=lig_number,
                                                      enumeration=0,
                                                      molecule=mol,
                                                      mol_type=self._LP.TYPE_RDKIT))

        list_ligands_enum = []
        for lig_number, mol in enumerate(Chem.SDMolSupplier(self.ligand_enum_sdf, removeHs=False)):
            list_ligands_enum.append(Ligand(smile=to_smiles(mol),
                                            original_smile=to_smiles(mol),
                                            ligand_number=lig_number + 50,   # + 50 to avoid lig_number overlap
                                            enumeration=0,
                                            molecule=mol,
                                            mol_type=self._LP.TYPE_RDKIT))

        self.ligands_with_hydrogens = list_ligands_with_hydrogens
        self.ligands_enum = list_ligands_enum

    def test_AutoDockVina_ligand_pdbqt_generation(self):
        docker = AutodockVina(
            input_pools=["RDkit"],
            parameters=AutodockVinaParameters(
                prefix_execution="module load AutoDock_Vina",
                search_space=SearchSpace(
                    center_x=0,
                    center_y=0,
                    center_z=0,
                    size_x=0,
                    size_y=0,
                    size_z=0
                )
            )
        )

        docker._initialize_executors()

        # generate PDBQT files and store paths
        ligands_with_hydrogens_paths = []
        return_ligands_with_hydrogens_list = []
        for ligand in self.ligands_with_hydrogens:
            cur_path = (os.path.join(self._folder_dir, "".join(["ligand_with_hydrogens_", str(ligand.get_ligand_number()), ".pdbqt"])))
            return_ligands_with_hydrogens_list.append(docker._write_molecule_to_pdbqt(cur_path, ligand.get_molecule()))
            ligands_with_hydrogens_paths.append(cur_path)

        self.assertEqual(sum(return_ligands_with_hydrogens_list), 15)
        self.assertEqual(lines_in_file(ligands_with_hydrogens_paths[0]), 72)
        self.assertEqual(lines_in_file(ligands_with_hydrogens_paths[2]), 62)
        self.assertEqual(lines_in_file(ligands_with_hydrogens_paths[7]), 63)

        ligands_enum_paths = []
        return_ligands_enum_list = []
        for ligand in self.ligands_enum:
            cur_path = (os.path.join(self._folder_dir, "".join(["ligand_enum_", str(ligand.get_ligand_number()), ".pdbqt"])))
            return_ligands_enum_list.append(docker._write_molecule_to_pdbqt(cur_path, ligand.get_molecule()))
            ligands_enum_paths.append(cur_path)

        self.assertEqual(sum(return_ligands_enum_list), 21)
        self.assertEqual(lines_in_file(ligands_enum_paths[0]), 69)
        self.assertEqual(lines_in_file(ligands_enum_paths[9]), 63)
        self.assertEqual(lines_in_file(ligands_enum_paths[20]), 92)

    def test_AutoDockVina_docking(self):
        docker = AutodockVina(
            input_pools=["RDkit"],
            parameters=AutodockVinaParameters(
                parallelization=Parallelization(number_cores=1),
                number_poses=4,
                receptor_pdbqt_path=[self.receptor_path],
                seed=11,
                search_space=SearchSpace(
                    center_x=3.3,
                    center_y=11.5,
                    center_z=24.8,
                    size_x=15,
                    size_y=10,
                    size_z=10
                ),
                prefix_execution="module load AutoDock_Vina"
            )
        )
        # Add the ligands_with_hydrogens molecules
        docker.add_molecules(molecules=self.ligands_with_hydrogens[:4])
        # Add the ligands_enumeration molecules
        docker.add_molecules(molecules=self.ligands_enum[:4])
        docker.dock()

        # write out the result and check length
        out_path = attach_root_path("tests/junk/ADV_docked.sdf")
        docker.write_docked_ligands(path=out_path)
        self.assertEqual(lines_in_file(out_path), 2606)

        # test the docked ligand poses (make sure, the read-out is not affected by deletion, changes or additions)
        list_docked_ligands = docker.get_docked_ligands()
        del list_docked_ligands[0]
        list_docked_ligands = docker.get_docked_ligands()
        self.assertEqual(8, len(list_docked_ligands))

        # test conformers
        self.assertListEqual([conf.GetProp(self._ROE.SDF_TAG_SCORE) for conf in docker.get_docked_ligands()[0].get_conformers()],
                             ["-9.1", "-8.1", "-7.9", "-7.8"])

        # Each row represents the top 4 poses of 1 ligand.
        # First 4 rows = ligands_with_hydrogens, last 4 rows = ligands_enumeration
        self.assertListEqual(docker.get_scores(best_only=False),
                             [-9.1, -8.1, -7.9, -7.8,
                              -9.0, -8.9, -8.1, -8.0,
                              -8.9, -8.2, -7.8, -7.5,
                              -9.2, -8.1, -7.9, -7.7,
                              -9.2, -8.9, -8.7, -8.2,
                              -9.2, -9.2, -8.4, -8.4,
                              -9.3, -8.4, -8.2,
                              -9.4, -9.0, -8.2, -8.0])
        self.assertListEqual(docker.get_scores(best_only=True),
                             [-9.1, -9.0, -8.9, -9.2, -9.2, -9.2, -9.3, -9.4])

        # test result parser (scores)
        path_scores_all = attach_root_path("tests/junk/ADV_result_parser_test_all.csv")
        docker.write_result(path=path_scores_all,
                            mode=self._CE.OUTPUT_MODE_ALL)
        self.assertEqual(lines_in_file(path_scores_all), 32)
        path_scores_best_per_enumeration = attach_root_path("tests/junk/ADV_result_parser_test_bpe.csv")
        docker.write_result(path=path_scores_best_per_enumeration,
                            mode=self._CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertEqual(lines_in_file(path_scores_best_per_enumeration), 9)

        # test result parser (poses)
        path_poses_best_per_enumeration = attach_root_path("tests/junk/ADV_result_parser_test_bpe.sdf")
        docker.write_docked_ligands(path=path_poses_best_per_enumeration,
                                    mode=self._CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertEqual(lines_in_file(path_poses_best_per_enumeration), 672)

    def test_AutoDockVina_docking_parallelized(self):
        docker = AutodockVina(
            input_pools=["RDkit"],
            parameters=AutodockVinaParameters(
                parallelization=Parallelization(number_cores=5),
                number_poses=4,
                receptor_pdbqt_path=[self.receptor_path],
                seed=9,
                search_space=SearchSpace(
                    center_x=3.3,
                    center_y=11.5,
                    center_z=24.8,
                    size_x=15,
                    size_y=10,
                    size_z=10
                ),
                prefix_execution="module load AutoDock_Vina"
            )
        )

        # Add the ligands_with_hydrogens molecules
        docker.add_molecules(molecules=self.ligands_with_hydrogens)

        # Add the ligands_enumeration molecules
        docker.add_molecules(molecules=self.ligands_enum)
        docker.dock()

        # write out the result and check length
        out_path = attach_root_path("tests/junk/ADV_docked_parallelized.sdf")
        docker.write_docked_ligands(path=out_path)
        self.assertEqual(lines_in_file(out_path), 8840)

        # test the docked ligand poses (make sure, the read-out is not affected by deletion, changes or additions)
        list_docked_ligands = docker.get_docked_ligands()
        del list_docked_ligands[0]
        list_docked_ligands = docker.get_docked_ligands()
        self.assertEqual(36, len(list_docked_ligands))

        # test conformers
        self.assertListEqual(
            [conf.GetProp(self._ROE.SDF_TAG_SCORE) for conf in docker.get_docked_ligands()[0].get_conformers()],
            ["-9.1", "-8.1", "-7.9", "-7.3"])
        self.assertListEqual(docker.get_scores(best_only=False),
                                [-9.1, -8.1, -7.9, -7.3,
                                 -9.0, -8.7, -8.1, -8.1,
                                 -8.9, -8.2, -7.8, -7.4,
                                 -9.2, -8.1, -7.9, -7.7,
                                 -9.3, -8.4, -8.0, -7.9,
                                 -8.9, -7.9, -7.9, -6.2,
                                 -9.7, -8.7, -8.4, -7.5,
                                 -9.3, -8.6, -8.3, -7.3,
                                 -10.0, -8.8, -7.6,
                                 -9.3, -8.7, -8.1, -6.8,
                                 -9.7, -9.6, -9.0, -7.8,
                                 -8.3, -7.8, -7.1,
                                 -9.1, -7.4, -6.8,
                                 -9.5, -9.1, -9.0, -8.4,
                                 -10.7, -9.6, -9.2, -9.1,
                                 -9.2, -8.9, -8.9, -8.5,
                                 -9.3, -9.3, -8.5, -8.5,
                                 -9.3, -8.4, -7.7, -7.3,
                                 -9.5, -8.3, -8.0, -7.8,
                                 -9.6, -9.2, -9.1, -8.9,
                                 -9.3, -8.2, -8.1, -6.9,
                                 -10.1, -9.0, -7.8, -7.8,
                                 -9.4, -8.5, -7.4, -6.6,
                                 -10.3, -9.1, -8.2,
                                 -9.4, -8.6, -8.5, -7.6,
                                 'NA', 'NA', 'NA',
                                 -9.5, -9.1, -9.0, -8.5,
                                 -9.8, -9.7, -9.5, -8.1,
                                 'NA', 'NA', 'NA',
                                 'NA', 'NA', 'NA'])
        self.assertListEqual(docker.get_scores(best_only=True),
                             [-9.1, -9.0, -8.9, -9.2, -9.3, -8.9, -9.7, -9.3,
                              -10.0, -9.3, -9.7, -8.3, -9.1, -9.5, -10.7, -9.2,
                              -9.3, -9.3, -9.5, -9.6, -9.3, -10.1, -9.4, -10.3,
                              -9.4, 'NA', 'NA', 'NA', -9.5, -9.8, 'NA', 'NA', 'NA',
                              'NA', 'NA', 'NA'])

        # test result parser (scores)
        path_scores_all = attach_root_path("tests/junk/ADV_parallelized_result_parser_test_all.csv")
        docker.write_result(path=path_scores_all,
                            mode=self._CE.OUTPUT_MODE_ALL)
        self.assertEqual(lines_in_file(path_scores_all), 105)
        path_scores_best_per_enumeration = attach_root_path("tests/junk/ADV_parallelized_result_parser_test_bpe.csv")
        docker.write_result(path=path_scores_best_per_enumeration,
                            mode=self._CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertEqual(lines_in_file(path_scores_best_per_enumeration), 28)

        # test result parser (poses)
        path_poses_best_per_enumeration = attach_root_path("tests/junk/ADV_parallelized_result_parser_test_bpe.sdf")
        docker.write_docked_ligands(path=path_poses_best_per_enumeration,
                                    mode=self._CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertEqual(lines_in_file(path_poses_best_per_enumeration), 2297)
