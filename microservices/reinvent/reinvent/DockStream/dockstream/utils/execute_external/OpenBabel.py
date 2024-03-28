import os
import sys
from dockstream.utils.enums.OpenBabel_enums import OpenBabelExecutablesEnum
from dockstream.utils.execute_external.execute import ExecutorBase

EE = OpenBabelExecutablesEnum()


class OpenBabelExecutor(ExecutorBase):
    """For the execution of the "obabel" binary."""

    def __init__(self):
        # in case the environment is not activated, add the path to the binary here
        obabel_location = os.path.dirname(sys.executable)
        super().__init__(prefix_execution=None, binary_location=obabel_location)

    def execute(self, command: str, arguments: list, check=True, location=None):
        # check, whether a proper executable is provided
        if command not in [EE.OBABEL]:
            raise ValueError("Parameter command must be an element of the internal OpenBabel executable list.")

        return super().execute(command=command,
                               arguments=arguments,
                               check=check,
                               location=location)

    def is_available(self):
        # unfortunately, "obabel" does not return a meaningful return value (always '1'), so instead try to parse
        # the "stdout" of the standard message; note, that "OpenBabel" is part of the environment and should always work
        try:
            result = self.execute(command=EE.OBABEL,
                                  arguments=[],
                                  check=False)
            if EE.OBABEL_IDENTIFICATION_STRING in result.stdout:
                return True
            return False
        except Exception as e:
            return False
