import unittest
import os
import shutil
from rdkit import Chem

from dockstream.containers.target_preparation_container import TargetPreparationContainer
from dockstream.core.AutodockVina.AutodockVina_target_preparator import AutodockVinaTargetPreparator

from dockstream.utils.enums.AutodockVina_enums import AutodockTargetPreparationEnum
from dockstream.utils.enums.OpenBabel_enums import OpenBabelExecutablesEnum
from dockstream.utils.execute_external.OpenBabel import OpenBabelExecutor

from tests.tests_paths import PATHS_1UYD, PATH_AUTODOCKVINA_EXAMPLES
from dockstream.utils.files_paths import attach_root_path


class Test_AutoDockVina_target_preparation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._TE = AutodockTargetPreparationEnum()
        cls._EE = OpenBabelExecutablesEnum()
        cls._original_path = os.getcwd()
        cls._working_dir = os.path.join(attach_root_path(PATH_AUTODOCKVINA_EXAMPLES.TARGET_PREP_FOLDER))

        cls._OpenBabel_executor = OpenBabelExecutor()
        if not cls._OpenBabel_executor.is_available():
            raise Exception("Could not load OpenBabel backend - abort.")

        # generate temporary folder and set current path there
        if not os.path.isdir(cls._working_dir):
            os.makedirs(cls._working_dir)
        os.chdir(cls._working_dir)

    def setUp(self) -> None:
        self.target = Chem.MolFromPDBFile(attach_root_path(PATHS_1UYD.TARGET_APO_PDB))
        self._conf = {self._TE.TARGETPREP: {
            self._TE.INPUT_PATH: attach_root_path(PATHS_1UYD.TARGET_APO_PDB),
            self._TE.FIX: {self._TE.FIX_ENABLED: False},
            self._TE.RUNS: [
                {
                    self._TE.RUNS_BACKEND: self._TE.RUNS_BACKEND_AUTODOCKVINA,
                    self._TE.RUNS_OUTPUT: {
                        self._TE.RECEPTOR_PATH: attach_root_path(PATH_AUTODOCKVINA_EXAMPLES.RECEPTOR)
                    },
                    self._TE.RUNS_PARAM: {
                        self._TE.PH: 7.4
                    }
                }
            ]
        }}

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls._original_path)
        if os.path.isdir(cls._working_dir):
            shutil.rmtree(cls._working_dir)

    def test_pdb2pdbqt(self):
        # specify PDBQT path
        pdbqt_path = os.path.join(self._working_dir, "test_pdbqt_gen.pdbqt")
        self._conf[self._TE.TARGETPREP][self._TE.RUNS][0][self._TE.RUNS_OUTPUT][self._TE.RECEPTOR_PATH] = pdbqt_path

        # write it out using openbabel
        conf = TargetPreparationContainer(conf=self._conf)
        prep = AutodockVinaTargetPreparator(conf=conf, target=self.target)
        prep._export_as_pdb2pdbqt(path=pdbqt_path)

        # check file size in bytes
        self.assertAlmostEqual(282410, os.stat(pdbqt_path).st_size, places=1)

    def test_ref_compound_box(self):
        self._conf[self._TE.TARGETPREP][self._TE.RUNS][0][self._TE.RUNS_PARAM][self._TE.EXTRACT_BOX] = {
            self._TE.EXTRACT_BOX_REFERENCE_LIGAND_PATH: attach_root_path(PATHS_1UYD.LIGAND_PU8_SDF),
            self._TE.EXTRACT_BOX_REFERENCE_LIGAND_FORMAT: self._TE.EXTRACT_BOX_REFERENCE_LIGAND_FORMAT_SDF
        }

        # write it out using openbabel
        conf = TargetPreparationContainer(conf=self._conf)
        prep = AutodockVinaTargetPreparator(conf=conf, target=self.target)
        x_coords, y_coords, z_coords = prep._extract_box()

        self.assertEqual(len(x_coords), 28)
        self.assertListEqual([4.403, 5.122, 5.091], x_coords[:3])
        self.assertListEqual([15.528, 15.084, 13.786], y_coords[:3])
        self.assertListEqual([26.579, 25.453, 24.846], z_coords[:3])
