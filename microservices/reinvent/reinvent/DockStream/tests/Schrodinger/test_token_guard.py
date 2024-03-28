import unittest

from dockstream.core.Schrodinger.license_token_guard import SchrodingerLicenseTokenGuard

from tests.tests_paths import PATH_SCHRODINGER_EXAMPLES
from dockstream.utils.files_paths import attach_root_path


class Test_Schrodinger_token_guard(unittest.TestCase):

    def test_token_guard_with_file(self):
        f = open(attach_root_path(PATH_SCHRODINGER_EXAMPLES.LICADMIN_FILE), 'r')
        lines = [line.rstrip("\n") for line in f.readlines()]
        f.close()

        # this works, as 144 tokens are available
        token_guard = SchrodingerLicenseTokenGuard(token_pools={"GLIDE_SP_DOCKING": 118})
        self.assertTrue(token_guard._check_licstat_output(lines))

        # this does not work, as only 144 tokens are available
        token_guard = SchrodingerLicenseTokenGuard(token_pools={"GLIDE_SP_DOCKING": 145})
        self.assertFalse(token_guard._check_licstat_output(lines))

        # this does not work, as only 200 tokens are available for GLIDE_XP_DOCKING
        token_guard = SchrodingerLicenseTokenGuard(token_pools={"GLIDE_SP_DOCKING": 118,
                                                                "GLIDE_XP_DOCKING": 201})
        self.assertFalse(token_guard._check_licstat_output(lines))

    def test_token_guard_with_execution(self):
        # this should work (until no more licenses are available at this point in time)
        token_guard = SchrodingerLicenseTokenGuard(
            prefix_execution="module load schrodinger/2019-4",
            token_pools={"GLIDE_SP_DOCKING": 8},
            wait_interval_seconds=10,
            wait_limit_seconds=60
        )
        self.assertTrue(token_guard.guard())

        # this will fail
        token_guard = SchrodingerLicenseTokenGuard(
            prefix_execution="module load schrodinger/2019-4",
            token_pools={"GLIDE_SP_DOCKING": 9000},
            wait_interval_seconds=3,
            wait_limit_seconds=10
        )
        self.assertFalse(token_guard.guard())
