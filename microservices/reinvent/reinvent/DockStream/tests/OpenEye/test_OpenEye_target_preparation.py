import unittest
import os
import shutil

import openeye.oedocking as oedocking

from dockstream.core.OpenEye.OpenEye_target_preparator import OpenEyeTargetPreparator
from dockstream.containers.target_preparation_container import TargetPreparationContainer

from dockstream.utils.enums.OpenEye_enums import OpenEyeTargetPreparationEnum
from tests.tests_paths import PATHS_1UYD, MAIN_CONFIG
from dockstream.utils.files_paths import attach_root_path


class Test_OpenEye_target_preparation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._original_path = os.getcwd()
        cls._working_dir = os.path.join(attach_root_path("tests/junk/OpenEye_target_prep"))
        if "OE_LICENSE" in MAIN_CONFIG:
            os.environ["OE_LICENSE"] = MAIN_CONFIG["OE_LICENSE"]

        # generate temporary folder and set current path there
        if not os.path.isdir(cls._working_dir):
            os.makedirs(cls._working_dir)
        os.chdir(cls._working_dir)

        cls._TE = OpenEyeTargetPreparationEnum()

    def setUp(self):
        self.target_path = attach_root_path(PATHS_1UYD.TARGET_APO_PDB)

        self._conf = {self._TE.TARGETPREP: {
            self._TE.INPUT_PATH: attach_root_path(PATHS_1UYD.TARGET_APO_PDB),
            self._TE.FIX: {self._TE.FIX_ENABLED: False},
            self._TE.RUNS: [
                {
                    self._TE.RUNS_BACKEND: self._TE.RUNS_BACKEND_OPENEYE,
                    self._TE.RUNS_OUTPUT: {self._TE.OUTPUT_RECEPTORPATH: ""},
                    self._TE.RUNS_PARAM: {},
                    self._TE.CAVITY: {
                        self._TE.CAVITY_METHOD: self._TE.CAVITY_METHOD_REFERENCE,
                        self._TE.CAVITY_REFERENCE_PATH: attach_root_path(PATHS_1UYD.LIGAND_PU8_SDF),
                        self._TE.CAVITY_REFERENCE_FORMAT: self._TE.CAVITY_REFERENCE_FORMAT_SDF
                    }
                },
                {
                    self._TE.RUNS_BACKEND: self._TE.RUNS_BACKEND_OPENEYE,
                    self._TE.RUNS_OUTPUT: {self._TE.OUTPUT_RECEPTORPATH: ""},
                    self._TE.RUNS_PARAM: {},
                    self._TE.CAVITY: {
                        self._TE.CAVITY_METHOD: self._TE.CAVITY_METHOD_HINT,
                        self._TE.CAVITY_HINT_COORDINATES: [
                            "0.234", 2.4, -0.89
                        ]
                    }
                },
                {
                    self._TE.RUNS_BACKEND: self._TE.RUNS_BACKEND_OPENEYE,
                    self._TE.RUNS_OUTPUT: {self._TE.OUTPUT_RECEPTORPATH: ""},
                    self._TE.RUNS_PARAM: {},
                    self._TE.CAVITY: {
                        self._TE.CAVITY_METHOD: self._TE.CAVITY_METHOD_BOX,
                        self._TE.CAVITY_BOX_LIMITS: [
                            '22', 15, 9, 4, 5, 3
                        ]
                    }
                }
            ]
        }}

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls._original_path)
        if os.path.isdir(cls._working_dir):
            shutil.rmtree(cls._working_dir)

    @staticmethod
    def _outer_contour_volume(receptor):
        negImagePot = oedocking.OEReceptorGetNegativeImageGrid(receptor)
        outerContourLevel = oedocking.OEReceptorGetOuterContourLevel(receptor)
        outerCount = 0
        for i in range(negImagePot.GetSize()):
            if negImagePot[i] >= outerContourLevel:
                outerCount += 1
        countToVolume = pow(negImagePot.GetSpacing(), 3)
        return outerCount * countToVolume

    def test_initialization(self):
        conf = TargetPreparationContainer(conf=self._conf)

        # initialize with a PDB file path ("reference" method)
        prep_reference = OpenEyeTargetPreparator(conf=conf, target=self.target_path, run_number=0)
        prep_reference.specify_cavity()
        reference = prep_reference.get_target()
        self.assertTrue(oedocking.OEReceptorHasBoundLigand(reference))
        self.assertEqual(962.2963823322922, self._outer_contour_volume(reference))

        # initialize with a PDB file path ("hint" method)
        prep_hint = OpenEyeTargetPreparator(conf=conf, target=self.target_path, run_number=1)
        prep_hint.specify_cavity()
        hint = prep_hint.get_target()
        self.assertTrue(not oedocking.OEReceptorHasBoundLigand(hint))
        self.assertEqual(1784.1853447037763, self._outer_contour_volume(hint))

        # initialize with a PDB file path ("box" method)
        prep_box = OpenEyeTargetPreparator(conf=conf, target=self.target_path, run_number=2)
        prep_box.specify_cavity()
        box = prep_box.get_target()
        self.assertTrue(not oedocking.OEReceptorHasBoundLigand(box))
        self.assertEqual(593.7037567849528, self._outer_contour_volume(box))
