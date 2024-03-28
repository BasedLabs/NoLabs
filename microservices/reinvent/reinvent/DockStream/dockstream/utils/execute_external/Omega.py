from dockstream.utils.enums.Omega_enums import OmegaExecutablesEnum
from dockstream.utils.execute_external.execute import ExecutorBase

_OE = OmegaExecutablesEnum()


class OmegaExecutor(ExecutorBase):
    """For the execution of the "OMEGA"."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(prefix_execution=prefix_execution, binary_location=binary_location)

    def execute(self, command: str, arguments: list, check=True, location=None):
        # check, whether a proper executable is provided
        if command not in [_OE.OMEGA]:
            raise ValueError("Command must be a valid parameter in the internal OMEGA dictionary.")

        return super().execute(command=command,
                               arguments=arguments,
                               check=check,
                               location=location)

    def is_available(self):
        try:
            result = self.execute(command=_OE.OMEGA,
                                  arguments=[_OE.HELP_SIMPLE],
                                  check=False)

            if _OE.OMEGA_HELP_IDENTIFICATION_STRING in result.stdout:
                return True
            return False
        except Exception as e:
            return False