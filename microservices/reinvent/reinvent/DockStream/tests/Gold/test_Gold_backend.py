import unittest
import os
import rdkit.Chem as Chem
from tests.tests_paths import MAIN_CONFIG
if "CSDHOME" in MAIN_CONFIG:
    os.environ["CSDHOME"] = MAIN_CONFIG["CSDHOME"]

from dockstream.core.Gold.Gold_docker import Gold, GoldParameters, GoldFitnessFunction, GoldResponseValue
from dockstream.utils.enums.Gold_enums import GoldDockingConfigurationEnum, GoldLigandPreparationEnum

from tests.tests_paths import PATHS_1UYD, PATH_GOLD_EXAMPLES
from dockstream.utils.files_paths import attach_root_path
from dockstream.core.ligand.ligand import Ligand
from dockstream.utils.smiles import to_smiles

from dockstream.utils.files_paths import lines_in_file


class Test_Gold_backend(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._CE = GoldDockingConfigurationEnum()
        cls._LP = GoldLigandPreparationEnum()

        # specify absolute paths to the input file(s)
        cls.target_path = attach_root_path(PATH_GOLD_EXAMPLES.TARGETFILE)
        cls.ligands_sdf = attach_root_path(PATHS_1UYD.LIGANDS_WITH_ENUMERATION_SDF)

    def setUp(self):
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
        self.ligands = list_ligands

    def test_Gold_docking_fitness(self):
        docker = Gold(
            input_pools=["Corina_pool"],
            parameters=GoldParameters(
                prefix_execution="module load ccdc/2020.3.0",
                receptor_paths=[self.target_path],
                fitness_function=GoldFitnessFunction.PLP,
                diverse_solutions=(False, None, None),
                early_termination=True,
                autoscale=100.0,
                ndocks=3,
            )
        )
        docker.add_molecules(molecules=self.ligands)
        docker.dock()
        result = docker.get_result()

        # test dataframe output
        self.assertEqual(63, result.shape[0])
        self.assertEqual(7, result.shape[1])
        self.assertTrue(all(float(x) > 50 for x in list(result.iloc[::3, :]["score"])))

        # test ligand retrieval
        list_docked_ligands = docker.get_docked_ligands()
        del list_docked_ligands[0]
        list_docked_ligands = docker.get_docked_ligands()
        self.assertEqual(21, len(list_docked_ligands))

        # test score retrieval
        self.assertEqual(17, len(docker.get_scores(best_only=True)))
        all_scores = docker.get_scores(best_only=False)
        self.assertEqual(63, len(all_scores))
        self.assertTrue(all_scores[0] > all_scores[1] > all_scores[2])
        self.assertTrue(all_scores[3] > all_scores[4] > all_scores[5])
        self.assertTrue(all_scores[6] > all_scores[7] > all_scores[8])

    def test_Gold_docking_value(self):
        docker = Gold(
            input_pools=["Corina_pool"],
            parameters=GoldParameters(
                prefix_execution="module load ccdc/2020.3.0",
                receptor_paths=[self.target_path],
                fitness_function=GoldFitnessFunction.PLP,
                diverse_solutions=(False, None, None),
                early_termination=True,
                autoscale=100.0,
                ndocks=3,
                response_value=GoldResponseValue.VALUE
            )
        )
        docker.add_molecules(molecules=self.ligands)
        docker.dock()
        result = docker.get_result()

        # test dataframe output
        self.assertEqual(63, result.shape[0])
        self.assertEqual(7, result.shape[1])
        self.assertTrue(all(float(x) < -40 for x in list(result.iloc[::3, :]["score"])))

        # test ligand retrieval
        list_docked_ligands = docker.get_docked_ligands()
        del list_docked_ligands[0]
        list_docked_ligands = docker.get_docked_ligands()
        self.assertEqual(21, len(list_docked_ligands))

        # write out poses and check length
        out_path = attach_root_path("tests/junk/Gold_backend_docked_all.sdf")
        docker.write_docked_ligands(path=out_path, mode=self._CE.OUTPUT_MODE_ALL)
        self.assertGreater(lines_in_file(out_path), 13000)
        out_path = attach_root_path("tests/junk/Gold_backend_docked_best_per_enumeration.sdf")
        docker.write_docked_ligands(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertGreater(lines_in_file(out_path), 4000)
        out_path = attach_root_path("tests/junk/Gold_backend_docked_best_per_ligand.sdf")
        docker.write_docked_ligands(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERLIGAND)
        self.assertGreater(lines_in_file(out_path), 3000)

        # write out dataframe and check length
        out_path = attach_root_path("tests/junk/Gold_backend_docked_all.csv")
        docker.write_result(path=out_path, mode=self._CE.OUTPUT_MODE_ALL)
        self.assertEqual(lines_in_file(out_path), 64)
        out_path = attach_root_path("tests/junk/Gold_backend_docked_best_per_enumeration.csv")
        docker.write_result(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertEqual(lines_in_file(out_path), 22)
        out_path = attach_root_path("tests/junk/Gold_backend_docked_best_per_ligand.csv")
        docker.write_result(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERLIGAND)
        self.assertEqual(lines_in_file(out_path), 18)
