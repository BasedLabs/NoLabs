import unittest
import os
import openeye.oechem as oechem
import openeye.oeomega as oeomega

from dockstream.core.OpenEye.OpenEye_docker import OpenEye, OpenEyeParameters
from dockstream.core.Schrodinger.Glide_docker import Parallelization
from dockstream.utils.enums.OpenEye_enums import OpenEyeDockingConfigurationEnum, OpenEyeLigandPreparationEnum

from dockstream.utils.enums.target_preparation_enum import TargetPreparationEnum
from dockstream.utils.enums.OpenEye_enums import OpenEyeResultKeywordsEnum

from tests.tests_paths import PATHS_1UYD, PATH_OPENEYE_EXAMPLES, MAIN_CONFIG
from dockstream.utils.files_paths import attach_root_path

from dockstream.core.ligand.ligand import Ligand
from dockstream.utils.smiles import to_smiles
from dockstream.utils.translations.translation import OpenEyeMolToRDkitMol

from dockstream.utils.files_paths import lines_in_file


class Test_OpenEye_backend(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._CE = OpenEyeDockingConfigurationEnum()
        cls._LP = OpenEyeLigandPreparationEnum()
        cls._TE = TargetPreparationEnum()
        cls._OE = OpenEyeResultKeywordsEnum()
        cls._RK = OpenEyeResultKeywordsEnum()

        # specify absolute paths to the input file(s)
        cls.target_path = attach_root_path(PATH_OPENEYE_EXAMPLES.RECEPTOR)
        cls.ligands_sdf = attach_root_path(PATHS_1UYD.LIGANDS_SDF_WITH_HYDROGENS)
        cls.ligands_enumerated = attach_root_path(PATHS_1UYD.LIGANDS_WITH_ENUMERATION_SDF)

        if "OE_LICENSE" in MAIN_CONFIG:
            os.environ["OE_LICENSE"] = MAIN_CONFIG["OE_LICENSE"]

    def setUp(self):
        # load the ligands
        ifs = oechem.oemolistream()
        ifs.SetFormat(oechem.OEFormat_SDF)
        ifs.open(self.ligands_sdf)
        list_ligands = []
        for lig_number, mol in enumerate(ifs.GetOEMols()):
            list_ligands.append(Ligand(smile=to_smiles(OpenEyeMolToRDkitMol(oechem.OEMol(mol), bySMILES=False)),
                                       original_smile=to_smiles(OpenEyeMolToRDkitMol(oechem.OEMol(mol), bySMILES=False)),
                                       ligand_number=lig_number,
                                       enumeration=0,
                                       molecule=oechem.OEMol(mol),
                                       mol_type=self._LP.TYPE_OPENEYE))
        ifs.close()
        self.ligands = list_ligands

        # enumerated
        ifs = oechem.oemolistream()
        ifs.SetFormat(oechem.OEFormat_SDF)
        ifs.open(self.ligands_enumerated)
        list_ligands = []
        for mol in ifs.GetOEMols():
            name = mol.GetTitle().split(':')
            lig_number = int(name[0])
            enumeration = int(name[1])
            list_ligands.append(Ligand(smile=to_smiles(OpenEyeMolToRDkitMol(oechem.OEMol(mol), bySMILES=False)),
                                       original_smile=to_smiles(OpenEyeMolToRDkitMol(oechem.OEMol(mol), bySMILES=False)),
                                       ligand_number=lig_number,
                                       enumeration=enumeration,
                                       molecule=oechem.OEMol(mol),
                                       mol_type=self._LP.TYPE_OPENEYE))
        ifs.close()
        self.enumerated_ligands = list_ligands

    def test_OpenEye_docking(self):
        docker = OpenEye(
            input_pools=["OpenEye_pool"],
            parameters=OpenEyeParameters(
                receptor_paths=[self.target_path],
                scoring=self._CE.SCORING_CHEMGAUSS3,
                resolution=self._CE.RESOLUTION_LOW,
                number_poses=3))

        docker.add_molecules(molecules=self.ligands)
        docker.dock()
        result = docker.get_result()

        # test dataframe output
        self.assertEqual(45, result.shape[0])
        self.assertEqual(7, result.shape[1])
        self.assertListEqual([-44.038658142089844, -53.41557693481445, -49.43500518798828,
                              -55.462684631347656, -31.567951202392578, -62.65696716308594,
                              -50.2397346496582, -45.879608154296875, -52.403839111328125,
                              -43.59895324707031, -37.22993850708008, -42.212188720703125,
                              -40.391151428222656, -60.124046325683594, -48.219669342041016],
                             list(result.iloc[::3, :]["score"]))

        # test ligand retrieval
        list_docked_ligands = docker.get_docked_ligands()
        del list_docked_ligands[0]
        list_docked_ligands = docker.get_docked_ligands()
        self.assertEqual(15, len(list_docked_ligands))
        self.assertEqual(3, len(list_docked_ligands[0].get_conformers()))

        # test score retrieval
        self.assertListEqual([-44.038658142089844, -53.41557693481445, -49.43500518798828, -55.462684631347656,
                              -31.567951202392578, -62.65696716308594, -50.2397346496582, -45.879608154296875,
                              -52.403839111328125, -43.59895324707031, -37.22993850708008, -42.212188720703125,
                              -40.391151428222656, -60.124046325683594, -48.219669342041016],
                             docker.get_scores(best_only=True))
        self.assertListEqual([-44.038658142089844, -38.612022399902344, -36.850990295410156, -53.41557693481445,
                              -46.19596862792969, -44.27885437011719, -49.43500518798828, -46.197410583496094,
                              -45.514591217041016, -55.462684631347656, -51.72509765625, -48.43769073486328,
                              -31.567951202392578, -30.920413970947266, -29.198612213134766, -62.65696716308594,
                              -51.419281005859375, -49.86194610595703, -50.2397346496582, -49.48503112792969,
                              -48.48686981201172, -45.879608154296875, -43.21759033203125, -41.363975524902344,
                              -52.403839111328125, -51.7275390625, -49.210018157958984, -43.59895324707031,
                              -42.1680908203125, -41.48081970214844, -37.22993850708008, -36.744483947753906,
                              -36.66645050048828, -42.212188720703125, -42.16796875, -42.12660217285156,
                              -40.391151428222656, -39.34788131713867, -38.36418151855469, -60.124046325683594,
                              -60.07300567626953, -59.94157028198242, -48.219669342041016, -43.59868240356445,
                              -39.95855712890625],
                             docker.get_scores(best_only=False))

    def test_OpenEye_docking_parallelized(self):
        docker = OpenEye(
            input_pools=["OpenEye_pool"],
            parameters=OpenEyeParameters(
                receptor_paths=[self.target_path],
                scoring=self._CE.SCORING_CHEMGAUSS3,
                parallelization=Parallelization(number_cores=4),
                resolution=self._CE.RESOLUTION_LOW,
                number_poses=3))
        docker.add_molecules(molecules=self.enumerated_ligands)
        docker.dock()
        result = docker.get_result()

        # test dataframe output
        self.assertEqual(51, result.shape[0])
        self.assertEqual(7, result.shape[1])
        self.assertListEqual([-44.038658142089844, -53.41557693481445, -49.43500518798828, -55.462684631347656,
                              -31.567951202392578, -62.65696716308594, -50.2397346496582, -45.879608154296875,
                              -52.403839111328125, -43.59895324707031, -41.7513542175293, -44.71626281738281,
                              -42.75702667236328, -55.63909912109375, -48.219669342041016, -36.374725341796875,
                              -53.31231689453125],
                             list(result.iloc[::3, :]["score"]))

        # test ligand retrieval
        list_docked_ligands = docker.get_docked_ligands()
        del list_docked_ligands[0]
        list_docked_ligands = docker.get_docked_ligands()
        self.assertEqual(21, len(list_docked_ligands))
        self.assertEqual(3, len(list_docked_ligands[0].get_conformers()))

        # test score retrieval
        self.assertListEqual([-44.038658142089844, -53.41557693481445, -49.43500518798828, -55.462684631347656,
                              -31.567951202392578, -62.65696716308594, -50.2397346496582, -45.879608154296875,
                              -52.403839111328125, -43.59895324707031, -41.7513542175293, -44.71626281738281,
                              -42.75702667236328, -55.63909912109375, -48.219669342041016, -53.31231689453125, 'NA'],
                             docker.get_scores(best_only=True))
        self.assertListEqual([-44.038658142089844, -38.612022399902344, -36.850990295410156, -53.41557693481445,
                              -46.19596862792969, -44.27885437011719, -49.43500518798828, -46.197410583496094,
                              -45.514591217041016, -55.462684631347656, -51.72509765625, -48.43769073486328,
                              -31.567951202392578, -30.920413970947266, -29.198612213134766, -62.65696716308594,
                              -51.419281005859375, -49.86194610595703, -50.2397346496582, -49.48503112792969,
                              -48.48686981201172, -45.879608154296875, -43.21759033203125, -41.363975524902344,
                              -52.403839111328125, -51.7275390625, -49.210018157958984, -43.59895324707031,
                              -42.1680908203125, -41.48081970214844, -41.7513542175293, -40.39088439941406,
                              -40.03367233276367, -44.71626281738281, -44.36548614501953, -40.34904479980469,
                              -42.75702667236328, -42.50238800048828, -41.681358337402344, -55.63909912109375,
                              -54.506778717041016, -54.029319763183594, -48.219669342041016, -43.59868240356445,
                              -39.95855712890625, -36.374725341796875, -35.98877716064453, -35.98469924926758,
                              -53.31231689453125, -52.625831604003906, -49.69447326660156, 'NA'],
                             docker.get_scores(best_only=False))

        # write out poses and check length
        out_path = attach_root_path("tests/junk/OpenEye_backend_docked_all.sdf")
        docker.write_docked_ligands(path=out_path, mode=self._CE.OUTPUT_MODE_ALL)
        self.assertEqual(lines_in_file(out_path), 5637)
        out_path = attach_root_path("tests/junk/OpenEye_backend_docked_best_per_enumeration.sdf")
        docker.write_docked_ligands(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertEqual(lines_in_file(out_path), 1879)
        out_path = attach_root_path("tests/junk/OpenEye_backend_docked_best_per_ligand_.sdf")
        docker.write_docked_ligands(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERLIGAND)
        self.assertEqual(lines_in_file(out_path), 1799)

        # write out dataframe and check length
        out_path = attach_root_path("tests/junk/OpenEye_backend_docked_all.csv")
        docker.write_result(path=out_path, mode=self._CE.OUTPUT_MODE_ALL)
        self.assertEqual(lines_in_file(out_path), 52)
        out_path = attach_root_path("tests/junk/OpenEye_backend_docked_best_per_enumeration.csv")
        docker.write_result(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertEqual(lines_in_file(out_path), 18)
        out_path = attach_root_path("tests/junk/OpenEye_backend_docked_best_per_ligand.csv")
        docker.write_result(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERLIGAND)
        self.assertEqual(lines_in_file(out_path), 17)

    def test_failed_docking(self):
        # replace one of the molecules by a failing one
        builder = oeomega.OEConformerBuilder()
        invalid_smile = "C#CCCCn1c(Cc2cc(OC)c(OC)c(OC)c2Cl)nc2"
        ligand_molecule = oechem.OEMol()
        oechem.OESmilesToMol(ligand_molecule, invalid_smile)
        self.ligands[3].set_smile(invalid_smile)
        self.ligands[3].set_molecule(ligand_molecule)

        docker = OpenEye(
            input_pools=["OpenEye_pool"],
            parameters=OpenEyeParameters(
                receptor_paths=[self.target_path],
                scoring=self._CE.SCORING_CHEMGAUSS4,
                resolution=self._CE.RESOLUTION_HIGH,
                number_poses=1))
        docker.add_molecules(molecules=self.ligands)
        docker.dock()
        result = docker.get_result()

        # test dataframe and score output
        self.assertEqual(14, result.shape[0])
        self.assertEqual(7, result.shape[1])
        self.assertListEqual([-6.6772565841674805, -9.564581871032715, -9.748031616210938, -5.963813781738281,
                              -9.797786712646484, -9.564628601074219, -7.988088130950928, -8.37636661529541,
                              -7.984253406524658, -7.760062217712402, -8.318611145019531, -8.598983764648438,
                              -10.127073287963867, -7.587806701660156],
                             list(result["score"]))
        self.assertListEqual(docker.get_scores(best_only=True),
                             [-6.6772565841674805, -9.564581871032715, -9.748031616210938, self._RK.FIXED_VALUE_NA,
                              -5.963813781738281, -9.797786712646484, -9.564628601074219, -7.988088130950928,
                              -8.37636661529541, -7.984253406524658, -7.760062217712402, -8.318611145019531,
                              -8.598983764648438, -10.127073287963867, -7.587806701660156])

        # write out the result and check length
        out_path = attach_root_path("tests/junk/OpenEye_backend_failed_docking_docked.sdf")
        docker.write_docked_ligands(path=out_path)
        self.assertEqual(lines_in_file(out_path), 1605)
