import os
import unittest
from shutil import copyfile

import rdkit.Chem as Chem

from dockstream.core.Schrodinger.Glide_docker import Parallelization
from dockstream.core.rDock.rDock_docker import rDock, rDockParameters

from dockstream.utils.enums.rDock_enums import rDockDockingConfigurationEnum, rDockResultKeywordsEnum, \
    rDockRbdockOutputEnum, rDockTargetPreparationEnum
from dockstream.utils.enums.RDkit_enums import RDkitLigandPreparationEnum

from tests.tests_paths import PATHS_1UYD, PATH_RDOCK_EXAMPLES
from dockstream.utils.files_paths import attach_root_path
from dockstream.core.ligand.ligand import Ligand
from dockstream.utils.smiles import to_smiles

from dockstream.utils.files_paths import lines_in_file


# TODO: move "prefix_execution" to "config/tests_config/config.json"
class Test_rDock_backend(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._CE = rDockDockingConfigurationEnum()
        cls._LP = RDkitLigandPreparationEnum()
        cls._RK = rDockResultKeywordsEnum()
        cls._ROE = rDockRbdockOutputEnum()
        cls._TP = rDockTargetPreparationEnum()

        # specify absolute paths to the various input files
        cls.ligands_sdf = attach_root_path(PATHS_1UYD.LIGANDS_WITH_ENUMERATION_SDF)
        cls.ligands_sdf_aligned = attach_root_path(PATHS_1UYD.LIGANDS_SDF_ALIGNED)
        cls.ligands_sdf_tethered = attach_root_path(PATHS_1UYD.LIGANDS_SDF_ALIGNED_TETHERED)

        # update the PRM (WITHOUT using the target preparator, so that there is no unit test dependency)
        # need to make temporary folder
        cls._unit_tests_dir = attach_root_path("tests/junk/rDock_cavity")
        if not os.path.exists(cls._unit_tests_dir):
            os.makedirs(cls._unit_tests_dir)
        copyfile(attach_root_path(PATH_RDOCK_EXAMPLES.CAVITY_GRID),
                 os.path.join(cls._unit_tests_dir, "rbcavity_1UYD_updated_cav1.grd"))
        copyfile(attach_root_path(PATH_RDOCK_EXAMPLES.CAVITY_AS),
                 os.path.join(cls._unit_tests_dir, "rbcavity_1UYD_updated.as"))

        with open(attach_root_path(PATH_RDOCK_EXAMPLES.CAVITY_PRM), "r") as prm_file:
            prm_string = prm_file.read()
        prm_string = prm_string.replace(cls._TP.STRING_RECEPTOR_MOL2_PATH,
                                        attach_root_path(PATHS_1UYD.TARGET_APO_MOL2))
        prm_string = prm_string.replace(cls._TP.STRING_REFERENCE_LIGAND_SDF_PATH,
                                        attach_root_path(PATHS_1UYD.LIGAND_PU8_SDF))
        tmp_prm_file = os.path.join(cls._unit_tests_dir,
                                    "rbcavity_1UYD_updated.prm")
        with open(tmp_prm_file, "wt") as prm_file:
            prm_file.write(prm_string)
        cls.docking_prm = tmp_prm_file

    def test_rDock_free_docking(self):
        docker = rDock(
            input_pools=["RDkit"],
            parameters=rDockParameters(
                prefix_execution="module load rDock",
                number_poses=3,
                rbdock_prm_paths=[self.docking_prm]
            )
        )

        list_ligands = []
        for mol in Chem.SDMolSupplier(self.ligands_sdf, removeHs=False):
            name = mol.GetProp("_Name").split(':')
            lig_number = int(name[0])
            enumeration = int(name[1])
            list_ligands.append(Ligand(smile=to_smiles(mol),
                                       original_smile=to_smiles(mol),
                                       ligand_number=lig_number,
                                       enumeration=enumeration,
                                       molecule=mol,
                                       mol_type=self._LP.TYPE_RDKIT))
        docker.add_molecules(molecules=list_ligands)
        docker.dock()
        result = docker.get_result()

        # write out the result and check length
        out_path = attach_root_path("tests/junk/rDock_backend_free_docking_docked.sdf")
        docker.write_docked_ligands(path=out_path)
        self.assertEqual(lines_in_file(out_path), 12816)

        # test the docked ligand poses (make sure, the read-out is not affected by deletion, changes or additions)
        list_docked_ligands = docker.get_docked_ligands()
        del list_docked_ligands[0]
        list_docked_ligands = docker.get_docked_ligands()
        self.assertEqual(21, len(list_docked_ligands))

        # test ordering
        self.assertListEqual([conf.GetProp(self._ROE.SCORE) for conf in docker.get_docked_ligands()[0].get_conformers()],
                             ["-23.9037", "-23.0126", "-18.1864"])
        self.assertListEqual(docker.get_scores(best_only=False),
                             [-23.9037, -23.0126, -18.1864, -39.842, -29.4756, -9.15032, -38.5213, -33.3635, -26.8268,
                              -30.8283, -27.4429, -26.318, -29.9767, -23.96, -20.711, -41.6575, -40.5511, -39.2945,
                              -32.4125, -30.2108, -27.8381, -34.7737, -29.4112, -28.7756, -40.8851, -35.702, -30.3977,
                              -36.1445, -31.7061, -29.5818, -26.0086, -22.1258, -19.3283, -31.3891, -20.1427, 186.343,
                              -38.0461, -34.5426, -24.5158, -42.3642, -40.1041, -37.1018, -39.7718, -38.8633, -25.0144,
                              -44.6188, -36.3875, -35.7647, -60.3626, -57.7277, -50.6844, -17.3378, 61.6253, 151.511,
                              -25.3268, 10.583, 166.711, -17.2741, -9.43039, 62.6485, -23.1856, -1.28365, 18.5449])
        self.assertListEqual(docker.get_scores(best_only=True),
                             [-23.9037, -39.842, -38.5213, -30.8283, -29.9767, -41.6575, -32.4125, -34.7737, -40.8851,
                              -36.1445, -26.0086, -31.3891, -38.0461, -42.3642, -39.7718, -60.3626, -25.3268])

        # test result parser (scores)
        path_scores_all = attach_root_path("tests/junk/rDock_result_parser_test_all.csv")
        docker.write_result(path=path_scores_all,
                            mode=self._CE.OUTPUT_MODE_ALL)
        self.assertEqual(lines_in_file(path_scores_all), 64)
        path_scores_best_per_enumeration = attach_root_path("tests/junk/rDock_result_parser_test_best_per_enumeration.csv")
        docker.write_result(path=path_scores_best_per_enumeration,
                            mode=self._CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertEqual(lines_in_file(path_scores_best_per_enumeration), 22)
        path_scores_best_per_ligand = attach_root_path("tests/junk/rDock_result_parser_test_best_per_ligand.csv")
        docker.write_result(path=path_scores_best_per_ligand,
                            mode=self._CE.OUTPUT_MODE_BESTPERLIGAND)
        self.assertEqual(lines_in_file(path_scores_best_per_ligand), 18)

        # test result parser (poses)
        path_poses_all = attach_root_path("tests/junk/rDock_result_parser_test_all.sdf")
        docker.write_docked_ligands(path=path_poses_all,
                                    mode=self._CE.OUTPUT_MODE_ALL)
        self.assertEqual(lines_in_file(path_poses_all), 12816)
        path_poses_best_per_enumeration = attach_root_path("tests/junk/rDock_result_parser_test_best_per_enumeration.sdf")
        docker.write_docked_ligands(path=path_poses_best_per_enumeration,
                                    mode=self._CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertEqual(lines_in_file(path_poses_best_per_enumeration), 4272)
        path_poses_best_per_ligand = attach_root_path("tests/junk/rDock_result_parser_test_best_per_ligand.sdf")
        docker.write_docked_ligands(path=path_poses_best_per_ligand,
                                    mode=self._CE.OUTPUT_MODE_BESTPERLIGAND)
        self.assertEqual(lines_in_file(path_poses_best_per_ligand), 3378)

    def test_rDock_free_docking_parallelized(self):
        docker = rDock(
            input_pools=["RDkit"],
            parameters=rDockParameters(
                prefix_execution="module load rDock",
                parallelization=Parallelization(number_cores=4),
                number_poses=3,
                rbdock_prm_paths=[self.docking_prm]
            )
        )

        list_ligands = []
        for lig_number, mol in enumerate(Chem.SDMolSupplier(self.ligands_sdf)):
            list_ligands.append(Ligand(smile=to_smiles(mol),
                                       original_smile=to_smiles(mol),
                                       ligand_number=lig_number,
                                       enumeration=0,
                                       molecule=mol,
                                       mol_type=self._LP.TYPE_RDKIT))
        docker.add_molecules(molecules=list_ligands)
        docker.dock()

        # write out the result and check length
        out_path = attach_root_path("tests/junk/rDock_backend_free_docking_docked_parallelized.sdf")
        docker.write_docked_ligands(path=out_path)
        self.assertGreater(lines_in_file(out_path), 12200)

        # test the docked ligand poses (make sure, the read-out is not affected by deletion, changes or additions)
        list_docked_ligands = docker.get_docked_ligands()
        del list_docked_ligands[0]
        list_docked_ligands = docker.get_docked_ligands()
        self.assertEqual(21, len(list_docked_ligands))

        # test ordering (if not ordered, it would be ["-19.5908", "-16.0143", "-38.2726"])
        self.assertListEqual([conf.GetProp(self._ROE.SCORE) for conf in docker.get_docked_ligands()[0].get_conformers()],
                             ["-29.628", "-23.4947", "-18.6307"])

    def test_rDock_aligned_docking(self):
        docker = rDock(
            input_pools=["RDkit"],
            parameters=rDockParameters(
                prefix_execution="module load rDock",
                number_poses=1,
                rbdock_prm_paths=[self.docking_prm]
            )
        )
        list_ligands = []
        for lig_number, mol in enumerate(Chem.SDMolSupplier(self.ligands_sdf_aligned)):
            list_ligands.append(Ligand(smile=to_smiles(mol),
                                       original_smile=to_smiles(mol),
                                       ligand_number=lig_number,
                                       enumeration=0,
                                       molecule=mol,
                                       mol_type=self._LP.TYPE_RDKIT))
        docker.add_molecules(molecules=list_ligands)
        docker.dock()

        self.assertEqual(15, len(docker.get_docked_ligands()))
        self.assertListEqual([conf.GetProp(self._ROE.SCORE) for conf in docker.get_docked_ligands()[0].get_conformers()],
                             ["-21.2678"])

    def test_rDock_tethered_docking(self):
        docker = rDock(
            input_pools=["RDkit"],
            parameters=rDockParameters(
                prefix_execution="module load rDock",
                number_poses=1,
                rbdock_prm_paths=[self.docking_prm]
            )
        )
        list_ligands = []
        for lig_number, mol in enumerate(Chem.SDMolSupplier(self.ligands_sdf_tethered)):
            list_ligands.append(Ligand(smile=to_smiles(mol),
                                       original_smile=to_smiles(mol),
                                       ligand_number=lig_number,
                                       enumeration=0,
                                       molecule=mol,
                                       mol_type=self._LP.TYPE_RDKIT))
        docker.add_molecules(molecules=list_ligands)
        docker.dock()

        self.assertEqual(15, len(docker.get_docked_ligands()))
        self.assertListEqual([conf.GetProp(self._ROE.SCORE) for conf in docker.get_docked_ligands()[0].get_conformers()],
                             ["-20.2089"])
