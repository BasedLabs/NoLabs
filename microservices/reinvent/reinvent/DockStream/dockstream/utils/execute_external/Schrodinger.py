from dockstream.utils.enums.Schrodinger_enums import SchrodingerExecutablesEnum
from dockstream.utils.execute_external.execute import ExecutorBase


EE = SchrodingerExecutablesEnum()


class SchrodingerExecutor(ExecutorBase):
    """For the execution of "Schrodinger" binaries."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(prefix_execution=prefix_execution, binary_location=binary_location)

    def execute(self, command: str, arguments: list, check=True, location=None, tokenGuard=None):
        # if a token guard has been configured, check if enough license tokens are available for the execution
        if tokenGuard is not None:
            tokenGuard.guard()

        # for Schrodinger commands, we need to provide the path to the binaries; we expect an environment variable
        # called "SCHRODINGER" at this stage and replace the command with the appropriate call
        if command == EE.GLIDE:
            command = EE.GLIDE_CALL
        elif command == EE.SDCONVERT:
            command = EE.SDCONVERT_CALL
        elif command == EE.LIGPREP:
            command = EE.LIGPREP_CALL
        else:
            raise ValueError("Parameter command must be an dictionary of the internal Schrodinger executable list.")
        return super().execute(command=command,
                               arguments=arguments,
                               check=check,
                               location=location)

    def is_available(self):
        try:
            result = self.execute(command=EE.GLIDE,
                                  arguments=[EE.GLIDE_HELP],
                                  check=True)

            if EE.GLIDE_HELP_IDENTIFICATION_STRING in result.stdout:
                return True
            return False
        except Exception as e:
            return False
