from dockstream.utils.enums.AutodockVina_enums import AutodockVinaExecutablesEnum
from dockstream.utils.execute_external.execute import ExecutorBase


EE = AutodockVinaExecutablesEnum()


class AutodockVinaExecutor(ExecutorBase):

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(prefix_execution=prefix_execution, binary_location=binary_location)

    def execute(self, command: str, arguments: list, check=True, location=None):
        # check, whether a proper executable is provided
        if command not in [EE.VINA]:
            raise ValueError("Parameter command must be an dictionary of the internal AutoDock Vina executable list.")

        return super().execute(command=command,
                               arguments=arguments,
                               check=check,
                               location=None)

    def is_available(self):
        try:
            result = self.execute(command=EE.VINA,
                                  arguments=[EE.VINA_VERSION],
                                  check=True)

            # do not use return code, because this will not work when wrapped another time
            if EE.VINA_VERSION_IDENTIFICATION_STRING in result.stdout:
                return True
            return False
        except Exception as e:
            return False
