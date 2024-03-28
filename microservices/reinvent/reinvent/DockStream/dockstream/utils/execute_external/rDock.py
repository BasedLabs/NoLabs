import os
import subprocess
import shutil
from dockstream.utils.files_paths import move_up_directory
from dockstream.utils.enums.rDock_enums import rDockExecutablesEnum
from dockstream.utils.execute_external.execute import ExecutorBase


EE = rDockExecutablesEnum()


class rDockExecutor(ExecutorBase):
    """For the execution of "rDock" binaries."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(prefix_execution=prefix_execution, binary_location=binary_location)

    def get_root(self):
        if self._prefix_execution is not None:
            command = self._prefix_execution + " && which " + EE.RBDOCK
            result = subprocess.run([command],
                                    check=True,
                                    universal_newlines=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    shell=True)
            return move_up_directory(result.stdout, n=2)
        else:
            return move_up_directory(shutil.which(EE.RBDOCK), n=1)

    def set_env_vars(self):
        # get path to root and set it to "RBT_ROOT" (important for library imports)
        # also set "RBT_HOME" which is the home directory for the input files which HAVE TO BE specified
        # with relative paths
        os.environ[EE.RBT_ROOT] = self.get_root()
        os.environ[EE.RBT_HOME] = os.getcwd()

    def execute(self, command: str, arguments: list, check=True, location=None):
        # check, whether a proper executable is provided
        if command not in [EE.RBDOCK, EE.RBCAVITY]:
            raise ValueError("Parameter command must be an dictionary of the internal rDock executable list.")

        return super().execute(command=command,
                               arguments=arguments,
                               check=check,
                               location=None)

    def is_available(self):
        try:
            result = self.execute(command=EE.RBDOCK,
                                  arguments=[EE.RBDOCK_HELP],
                                  check=True)

            # do not use return code, because this will not work when wrapped another time
            if EE.RBDOCK_HELP_IDENTIFICATION_STRING in result.stdout:
                return True
            return False
        except Exception as e:
            return False
