from dockstream.utils.enums.Gold_enums import GoldExecutablesEnum
from dockstream.utils.execute_external.execute import ExecutorBase


EE = GoldExecutablesEnum()


class GoldExecutor(ExecutorBase):
    """For the execution of "Gold" binaries."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(prefix_execution=prefix_execution, binary_location=binary_location)

    def execute(self, command: str, arguments: list, check=True, location=None):
        # check, whether a proper executable is provided
        if command not in [EE.GOLD_AUTO]:
            raise ValueError("Parameter command must be an dictionary of the internal Gold executable list.")

        return super().execute(command=command,
                               arguments=arguments,
                               check=check,
                               location=location)

    def is_available(self):
        try:
            result = self.execute(command=EE.GOLD_AUTO,
                                  arguments=[EE.GOLD_AUTO_HELP],
                                  check=False)

            # do not use return code, because this will not work when wrapped another time
            if EE.GOLD_AUTO_HELP_IDENTIFICATION_STRING in result.stdout:
                return True
            return False
        except Exception as e:
            return False
