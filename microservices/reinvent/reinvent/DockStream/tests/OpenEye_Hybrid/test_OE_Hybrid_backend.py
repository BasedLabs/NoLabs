import unittest
import rdkit.Chem as Chem

from dockstream.core.OpenEyeHybrid.OpenEyeHybrid_docker import OpenEyeHybrid, OpenEyeHybridParameters
from dockstream.core.Schrodinger.Glide_docker import Parallelization
from dockstream.utils.enums.OE_Hybrid_enums import OpenEyeHybridDockingConfigurationEnum, \
                                               OpenEyeHybridLigandPreparationEnum, \
                                               OpenEyeHybridResultKeywordsEnum
from dockstream.utils.enums.OE_Hybrid_enums import OpenEyeHybridExecutablesEnum, OpenEyeHybridOutputKeywordsEnum
from dockstream.utils.enums.target_preparation_enum import TargetPreparationEnum

from tests.tests_paths import PATHS_1UYD, PATH_OPENEYEHYBRID_EXAMPLES
from dockstream.utils.files_paths import attach_root_path

from dockstream.core.ligand.ligand import Ligand
from dockstream.utils.smiles import to_smiles

from dockstream.utils.files_paths import lines_in_file


class Test_OpenEyeHybrid_backend(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._CE = OpenEyeHybridDockingConfigurationEnum()
        cls._LP = OpenEyeHybridLigandPreparationEnum()
        cls._TE = TargetPreparationEnum()
        cls._RK = OpenEyeHybridResultKeywordsEnum()
        cls._EE = OpenEyeHybridExecutablesEnum()
        cls._OE = OpenEyeHybridOutputKeywordsEnum()

        # specify absolute paths to the input file(s)
        cls.receptor_path = attach_root_path(PATH_OPENEYEHYBRID_EXAMPLES.RECEPTOR)
        cls.ligands_sdf = attach_root_path(PATHS_1UYD.LIGANDS_SDF_WITH_HYDROGENS)
        cls.ligands_enumerated = attach_root_path(PATHS_1UYD.LIGANDS_WITH_ENUMERATION_SDF)

    def setUp(self):
        # set up ligands with enumeration
        list_ligands = []
        for mol in Chem.SDMolSupplier(self.ligands_enumerated, removeHs=False):
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

    """ include the parameters in the DockingContainer run in the json"""
    def test_OpenEyeHybrid_docking(self):
        docker = OpenEyeHybrid(
            input_pools=["OpenEye"],
            parameters=OpenEyeHybridParameters(
                prefix_execution=self._EE.OE_HYBRID_MODULE_LOAD,
                receptor_paths=[self.receptor_path],
                resolution=self._CE.RESOLUTION_LOW,
                number_poses=3
            )
        )
        docker.add_molecules(molecules=self.ligands)
        docker.dock()
        result = docker.get_result()

        # test dataframe output
        self.assertEqual(51, result.shape[0])
        # test number of columns
        self.assertEqual(7, result.shape[1])
        # ensure number of ligands successfully docked is 17/21
        self.assertEqual(17, len(list(result.iloc[::3, :]["score"])))

        # ensure every 3rd score of the docked ligands is correct
        # the docking run was set to generate 3 poses per ligand so this should check 1 pose of each
        # successfully docked ligand
        self.assertEqual([-3.367321, -0.804976, -2.991787, -1.990613, -2.468151, -5.344979,
                          -1.844902, -2.727283, -5.000861, -2.41826, -3.696608, -4.858593,
                          -4.81919, -3.910051, -4.909737, -2.266009, -2.44015], list(result.iloc[::3, :]["score"]))

        # test ligand retrieval
        list_docked_ligands = docker.get_docked_ligands()
        del list_docked_ligands[0]
        self.assertEqual(20, len(list_docked_ligands))
        list_docked_ligands = docker.get_docked_ligands()
        self.assertEqual(21, len(list_docked_ligands))
        self.assertEqual(3, len(list_docked_ligands[0].get_conformers()))

        self.assertEqual([-3.367321, -0.804976, -2.991787, -1.990613, -2.468151, -5.344979, -1.844902,
                          -2.727283, -5.000861, -2.41826, -3.696608, -4.858593, -4.81919, -3.910051,
                          -4.909737, -2.44015, 'NA'], docker.get_scores(best_only=True))

        self.assertEqual([-3.367321, -2.674903, -2.359775,
                          -0.804976, -0.589188, -0.564974,
                          -2.991787, -2.735957, -2.33122,
                          -1.990613, -1.897298, -1.768375,
                          -2.468151, -2.105051, -0.788182,
                          -5.344979, -2.234832, -1.75003,
                          -1.844902, -1.606467, -1.501745,
                          -2.727283, -2.563828, -2.23263,
                          -5.000861, -3.259999, -2.687797,
                          -2.41826, -2.380332, -2.245145,
                          -3.696608, -3.592045, -2.570743,
                          -4.858593, -4.394607, -4.003824,
                          -4.81919, -4.264267, -3.589545,
                          -3.910051, -2.384474, -2.349108,
                          -4.909737, -4.847009, -4.824538,
                          -2.266009, -2.228133, -2.038738,
                          -2.44015, -2.382029, -2.309129,
                          'NA'], docker.get_scores(best_only=False))

        # write out poses and check length
        out_path = attach_root_path("tests/junk/OpenEyeHybrid_backend_docked_all.sdf")
        docker.write_docked_ligands(path=out_path, mode=self._CE.OUTPUT_MODE_ALL)
        self.assertEqual(lines_in_file(out_path), 6402)
        out_path = attach_root_path("tests/junk/OpenEyeHybrid_backend_docked_best_per_enumeration.sdf")
        docker.write_docked_ligands(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertEqual(lines_in_file(out_path), 2134)
        out_path = attach_root_path("tests/junk/OpenEyeHybrid_backend_docked_best_per_ligand.sdf")
        docker.write_docked_ligands(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERLIGAND)
        self.assertEqual(lines_in_file(out_path), 2039)

        # write out DataFrame and check length
        out_path = attach_root_path("tests/junk/OpenEyeHybrid_backend_docked_all.csv")
        docker.write_result(path=out_path, mode=self._CE.OUTPUT_MODE_ALL)
        self.assertEqual(lines_in_file(out_path), 52)
        out_path = attach_root_path("tests/junk/OpenEyeHybrid_backend_docked_best_per_enumeration.csv")
        docker.write_result(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertEqual(lines_in_file(out_path), 18)
        out_path = attach_root_path("tests/junk/OpenEyeHybrid_backend_docked_best_per_ligand.csv")
        docker.write_result(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERLIGAND)
        self.assertEqual(lines_in_file(out_path), 17)

    def test_OpenEyeHybrid_parallelized_docking(self):
        docker = OpenEyeHybrid(
            input_pools=["OpenEye"],
            parameters=OpenEyeHybridParameters(
                prefix_execution=self._EE.OE_HYBRID_MODULE_LOAD,
                parallelization=Parallelization(number_cores=4),
                receptor_paths=[self.receptor_path],
                resolution=self._CE.RESOLUTION_LOW,
                number_poses=3
            )
        )

        docker.add_molecules(molecules=self.ligands)
        docker.dock()

        result = docker.get_result()
        # test dataframe output
        self.assertEqual(51, result.shape[0])
        # test number of columns
        self.assertEqual(7, result.shape[1])
        # ensure number of ligands successfully docked is 17/21
        self.assertEqual(17, len(list(result.iloc[::3, :]["score"])))

        # ensure every 3rd score of the docked ligands is correct
        # the docking run was set to generate 3 poses per ligand so this should check 1 pose of each
        # successfully docked ligand
        self.assertEqual([-3.367321, -0.804976, -2.991787, -1.990613, -2.468151, -5.344979,
                          -1.844902, -2.727283, -5.000861, -2.41826, -3.696608, -4.858593,
                          -4.81919, -3.910051, -4.909737, -2.266009, -2.44015], list(result.iloc[::3, :]["score"]))

        # test ligand retrieval
        list_docked_ligands = docker.get_docked_ligands()
        del list_docked_ligands[0]
        self.assertEqual(20, len(list_docked_ligands))
        list_docked_ligands = docker.get_docked_ligands()
        self.assertEqual(21, len(list_docked_ligands))
        self.assertEqual(3, len(list_docked_ligands[0].get_conformers()))

        self.assertEqual([-3.367321, -0.804976, -2.991787, -1.990613, -2.468151, -5.344979, -1.844902,
                          -2.727283, -5.000861, -2.41826, -3.696608, -4.858593, -4.81919, -3.910051,
                          -4.909737, -2.44015, 'NA'], docker.get_scores(best_only=True))

        self.assertEqual([-3.367321, -2.674903, -2.359775,
                          -0.804976, -0.589188, -0.564974,
                          -2.991787, -2.735957, -2.33122,
                          -1.990613, -1.897298, -1.768375,
                          -2.468151, -2.105051, -0.788182,
                          -5.344979, -2.234832, -1.75003,
                          -1.844902, -1.606467, -1.501745,
                          -2.727283, -2.563828, -2.23263,
                          -5.000861, -3.259999, -2.687797,
                          -2.41826, -2.380332, -2.245145,
                          -3.696608, -3.592045, -2.570743,
                          -4.858593, -4.394607, -4.003824,
                          -4.81919, -4.264267, -3.589545,
                          -3.910051, -2.384474, -2.349108,
                          -4.909737, -4.847009, -4.824538,
                          -2.266009, -2.228133, -2.038738,
                          -2.44015, -2.382029, -2.309129,
                          'NA'], docker.get_scores(best_only=False))

        # write out poses and check length
        out_path = attach_root_path("tests/junk/OpenEyeHybrid_backend_parallelized_docked_all.sdf")
        docker.write_docked_ligands(path=out_path, mode=self._CE.OUTPUT_MODE_ALL)
        self.assertEqual(lines_in_file(out_path), 6402)
        out_path = attach_root_path("tests/junk/OpenEyeHybrid_backend_parallelized_docked_best_per_enumeration.sdf")
        docker.write_docked_ligands(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertEqual(lines_in_file(out_path), 2134)
        out_path = attach_root_path("tests/junk/OpenEyeHybrid_backend_parallelized_docked_best_per_ligand.sdf")
        docker.write_docked_ligands(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERLIGAND)
        self.assertEqual(lines_in_file(out_path), 2039)

        # write out DataFrame and check length
        out_path = attach_root_path("tests/junk/OpenEyeHybrid_backend_parallelized_docked_all.csv")
        docker.write_result(path=out_path, mode=self._CE.OUTPUT_MODE_ALL)
        self.assertEqual(lines_in_file(out_path), 52)
        out_path = attach_root_path("tests/junk/OpenEyeHybrid_backend_parallelized_docked_best_per_enumeration.csv")
        docker.write_result(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertEqual(lines_in_file(out_path), 18)
        out_path = attach_root_path("tests/junk/OpenEyeHybrid_backend_parallelized_docked_best_per_ligand.csv")
        docker.write_result(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERLIGAND)
        self.assertEqual(lines_in_file(out_path), 17)
