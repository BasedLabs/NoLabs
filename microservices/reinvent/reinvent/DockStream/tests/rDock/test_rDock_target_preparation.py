import unittest
import os
import shutil
from rdkit import Chem

from dockstream.containers.target_preparation_container import TargetPreparationContainer
from dockstream.core.rDock.rDock_target_preparator import rDockTargetPreparator

from dockstream.utils.enums.rDock_enums import rDockTargetPreparationEnum
from dockstream.utils.enums.rDock_enums import rDockExecutablesEnum, rDockResultKeywordsEnum
from dockstream.utils.execute_external.rDock import rDockExecutor

from tests.tests_paths import PATHS_1UYD, PATH_RDOCK_EXAMPLES
from dockstream.utils.files_paths import attach_root_path


class Test_rDock_target_preparation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._TE = rDockTargetPreparationEnum()
        cls._EE = rDockExecutablesEnum()
        cls._GK = rDockResultKeywordsEnum()
        cls._original_path = os.getcwd()
        cls._working_dir = os.path.join(attach_root_path("tests/junk/rDock_target_prep"))

        # TODO: add "prefix execution" setting to "config/tests_config/config.json"
        cls._rDock_executor = rDockExecutor(prefix_execution="module load rDock", binary_location=None)
        if not cls._rDock_executor.is_available():
            raise Exception("Could not load rDock backend - abort.")
        cls._rDock_executor.set_env_vars()

        # generate temporary folder and set current path there
        if not os.path.isdir(cls._working_dir):
            os.makedirs(cls._working_dir)
        os.chdir(cls._working_dir)

        # copy test PRM file to temporary folder
        path_cavity_prm = os.path.join(cls._working_dir, os.path.basename(PATH_RDOCK_EXAMPLES.CAVITY_PRM))
        path_ligand_sdf = attach_root_path(PATHS_1UYD.LIGAND_PU8_SDF)
        shutil.copyfile(attach_root_path(PATH_RDOCK_EXAMPLES.CAVITY_PRM), path_cavity_prm)

        cls._conf = {cls._TE.TARGETPREP: {
            cls._TE.INPUT_PATH: attach_root_path(PATHS_1UYD.TARGET_APO_PDB),
            cls._TE.FIX: {cls._TE.FIX_ENABLED: False},
            cls._TE.RUNS: [
                {
                    cls._TE.RUNS_BACKEND: cls._TE.RUNS_BACKEND_RDOCK,
                    cls._TE.RUNS_OUTPUT: {
                        cls._TE.RUNS_OUTPUT_DIRECTORY: cls._working_dir
                    },
                    cls._TE.RUNS_PARAM: {cls._TE.RUNS_PARAM_PREFIX_EXECUTION: "module load rDock"},
                    cls._TE.CAVITY: {
                        cls._TE.CAVITY_METHOD: cls._TE.CAVITY_METHOD_REFERENCE,
                        cls._TE.CAVITY_REFERENCE_PATH: path_ligand_sdf,
                        cls._TE.CAVITY_REFERENCE_FORMAT: cls._TE.CAVITY_REFERENCE_FORMAT_SDF,
                        cls._TE.CAVITY_PRMFILE: path_cavity_prm
                    }
                }
            ]
        }}

    def setUp(self):
        self.target = Chem.MolFromPDBFile(attach_root_path(PATHS_1UYD.TARGET_APO_PDB))
        self.target_mol2 = attach_root_path(PATHS_1UYD.TARGET_APO_MOL2)

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls._original_path)
        if os.path.isdir(cls._working_dir):
            shutil.rmtree(cls._working_dir)

    def test_mol2_generation(self):
        # specify Mol2 file path
        self._mol2_path = os.path.join(self._working_dir, "target.mol2")

        # write it out using openbabel (for atom type specification) in "private" method
        conf = TargetPreparationContainer(conf=self._conf)
        prep = rDockTargetPreparator(conf=conf, target=self.target)
        prep._export_as_mol2(path=self._mol2_path)

        # check file size in bytes
        self.assertAlmostEqual(411116/10000, os.stat(self._mol2_path).st_size/10000, places=1)

    def test_specify_cavity_external(self):
        conf = TargetPreparationContainer(conf=self._conf)
        prep = rDockTargetPreparator(conf=conf, target=self.target)
        result = prep.specify_cavity()

        self.assertEqual(result[self._GK.SPECIFYCAVITY_METADATA][self._GK.SPECIFYCAVITY_METADATA_TOTALVOLUME], 1025.75)
        self.assertEqual(result[self._GK.SPECIFYCAVITY_METADATA][self._GK.SPECIFYCAVITY_METADATA_SIZEINPOINTS], 8206)
        self.assertTrue(os.path.isfile(result[self._GK.SPECIFYCAVITY_BINARY_PATH]))
        self.assertAlmostEqual(1353664/10000,
                               os.stat(result[self._GK.SPECIFYCAVITY_BINARY_PATH]).st_size/10000,
                               places=1)
        self.assertTrue(os.path.isfile(result[self._GK.SPECIFYCAVITY_GRID_PATH]))

    def test_specify_cavity_internal(self):
        path_ligand_pdb = attach_root_path(PATHS_1UYD.LIGAND_PU8_PDB)
        conf = {self._TE.TARGETPREP: {
            self._TE.INPUT_PATH: attach_root_path(PATHS_1UYD.TARGET_APO_PDB),
            self._TE.FIX: {self._TE.FIX_ENABLED: False},
            self._TE.RUNS: [
                {
                    self._TE.RUNS_BACKEND: self._TE.RUNS_BACKEND_RDOCK,
                    self._TE.RUNS_OUTPUT: {
                        self._TE.RUNS_OUTPUT_DIRECTORY: self._working_dir
                    },
                    self._TE.RUNS_PARAM: {self._TE.RUNS_PARAM_PREFIX_EXECUTION: "module load rDock"},
                    self._TE.CAVITY: {
                        self._TE.CAVITY_METHOD: self._TE.CAVITY_METHOD_REFERENCE,
                        self._TE.CAVITY_REFERENCE_PATH: path_ligand_pdb,
                        self._TE.CAVITY_REFERENCE_FORMAT: self._TE.CAVITY_REFERENCE_FORMAT_PDB
                    }
                }
            ]
        }}
        conf = TargetPreparationContainer(conf=conf)
        prep = rDockTargetPreparator(conf=conf, target=self.target)
        result = prep.specify_cavity()

        self.assertEqual(result[self._GK.SPECIFYCAVITY_METADATA][self._GK.SPECIFYCAVITY_METADATA_TOTALVOLUME], 1025.75)
        self.assertEqual(result[self._GK.SPECIFYCAVITY_METADATA][self._GK.SPECIFYCAVITY_METADATA_SIZEINPOINTS], 8206)
        self.assertTrue(os.path.isfile(result[self._GK.SPECIFYCAVITY_BINARY_PATH]))
        self.assertAlmostEqual(1353664/10000,
                               os.stat(result[self._GK.SPECIFYCAVITY_BINARY_PATH]).st_size/10000,
                               places=1)
