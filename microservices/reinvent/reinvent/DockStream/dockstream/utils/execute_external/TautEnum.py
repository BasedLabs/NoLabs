from dockstream.utils.enums.taut_enum_enums import TautEnumEnum
from dockstream.utils.execute_external.execute import ExecutorBase


EE = TautEnumEnum()


class TautEnumExecutor(ExecutorBase):
    """For the execution of the "TautEnum" binary."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(prefix_execution=prefix_execution, binary_location=binary_location)

    def execute(self, command: str, arguments: list, check=True, location=None):
        # check, whether a proper executable is provided
        if command not in [EE.TAUTENUM]:
            raise ValueError("Parameter command must be an dictionary of the internal TautEnum executable list.")

        return super().execute(command=command,
                               arguments=arguments,
                               check=check,
                               location=location)

    def is_available(self):
        # unfortunately, "TautEnum" does not seem to return a meaningful return value, so instead try to parse
        # the "stdout" of the help message
        try:
            result = self.execute(command=EE.TAUTENUM,
                                  arguments=[EE.TAUTENUM_HELP],
                                  check=False)
            if EE.TAUTENUM_HELP_IDENTIFICATION_STRING in result.stdout:
                return True
            return False
        except Exception as e:
            return False
