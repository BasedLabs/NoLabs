import os
import shutil
import unittest

from tests.tests_paths import MAIN_CONFIG
if "CSDHOME" in MAIN_CONFIG:
    os.environ["CSDHOME"] = MAIN_CONFIG["CSDHOME"]

from dockstream.utils.enums.Gold_enums import GoldTargetKeywordEnum
from dockstream.containers.target_preparation_container import TargetPreparationContainer
from dockstream.core.Gold.Gold_target_preparator import GoldTargetPreparator

from tests.tests_paths import PATHS_1UYD
from dockstream.utils.files_paths import attach_root_path


class Test_Gold_target_preparation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._original_path = os.getcwd()
        cls._working_dir = os.path.join(attach_root_path("tests/junk/Gold_target_prep"))

        # generate temporary folder and set current path there
        if not os.path.isdir(cls._working_dir):
            os.makedirs(cls._working_dir)
        os.chdir(cls._working_dir)

        cls._TE = GoldTargetKeywordEnum()

    def setUp(self):
        self.target_path = attach_root_path(PATHS_1UYD.TARGET_APO_PDB)

        self._conf = {self._TE.TARGETPREP: {
            self._TE.INPUT_PATH: attach_root_path(PATHS_1UYD.TARGET_APO_PDB),
            self._TE.FIX: {self._TE.FIX_ENABLED: False},
            self._TE.RUNS: [
                {
                    self._TE.RUNS_BACKEND: self._TE.RUNS_BACKEND_GOLD,
                    self._TE.RUNS_OUTPUT: {},
                    self._TE.RUNS_PARAM: {},
                    self._TE.CAVITY: {
                        self._TE.CAVITY_METHOD: self._TE.CAVITY_METHOD_REFERENCE,
                        self._TE.CAVITY_REFERENCE_PATH: attach_root_path(PATHS_1UYD.LIGAND_PU8_PDB),
                        self._TE.CAVITY_REFERENCE_DISTANCE: 7.0
                    }
                }
            ]
        }}

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls._original_path)
        if os.path.isdir(cls._working_dir):
            shutil.rmtree(cls._working_dir)

    def test_initialization(self):
        conf = TargetPreparationContainer(conf=self._conf)

        # initialize with a PDB file path ("reference" method)
        prep_reference = GoldTargetPreparator(conf=conf, target=self.target_path, run_number=0)
        prep_reference.specify_cavity()
        prep_reference.write_target(path=os.path.join(self._working_dir, "Gold_binding_site.pkl"))
        dictionary = prep_reference._target_dict
        self.assertEqual(len(dictionary[self._TE.TARGET_PDB]), 1925)
        self.assertEqual(len(dictionary[self._TE.REFERENCE_LIGAND]), 28)
        self.assertEqual(dictionary[self._TE.CAVITY_REFERENCE_DISTANCE], 7.0)
