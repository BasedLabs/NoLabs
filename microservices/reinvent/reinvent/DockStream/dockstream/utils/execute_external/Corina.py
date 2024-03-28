from dockstream.utils.enums.Corina_enums import CorinaExecutablesEnum
from dockstream.utils.execute_external.execute import ExecutorBase


EE = CorinaExecutablesEnum()


class CorinaExecutor(ExecutorBase):
    """For the execution of the "Corina" binary."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(prefix_execution=prefix_execution, binary_location=binary_location)

    def execute(self, command: str, arguments: list, check=True, location=None):
        # check, whether a proper executable is provided
        if command not in [EE.CORINA]:
            raise ValueError("Parameter command must be an dictionary of the internal Corina executable list.")

        return super().execute(command=command,
                               arguments=arguments,
                               check=check,
                               location=location)

    def is_available(self):
        # unfortunately, "Corina" does not return a meaningful return value (always '1'), so instead try to parse
        # the "stderr" of the help message
        try:
            result = self.execute(command=EE.CORINA,
                                  arguments=[EE.CORINA_HELP],
                                  check=False)
            if EE.CORINA_HELP_IDENTIFICATION_STRING in result.stderr:
                return True
            return False
        except Exception as e:
            return False
