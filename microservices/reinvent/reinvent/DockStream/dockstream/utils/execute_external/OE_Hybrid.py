from dockstream.utils.enums.OE_Hybrid_enums import OpenEyeHybridExecutablesEnum
from dockstream.utils.execute_external.execute import ExecutorBase

EE = OpenEyeHybridExecutablesEnum()


class OpenEyeHybridExecutor(ExecutorBase):
    """For the execution of OpenEye Hybrid docking which optimizes enrichment"""

    def __init__(self, prefix_execution="module load oedocking", binary_location=None):
        super().__init__(prefix_execution=prefix_execution, binary_location=binary_location)

    def execute(self, command: str, arguments: list, check=True, location=None):
        # check, whether a proper executable is provided
        if command not in [EE.HYBRID]:
            raise ValueError("Parameter command must be a dictionary of the internal OpenEye Hybrid executable list.")

        return super().execute(command=command,
                               arguments=arguments,
                               check=check,
                               location=location)

    def is_available(self):
        try:
            result = self.execute(command=EE.HYBRID,
                                  arguments=[EE.OE_HYBRID_HELP_SIMPLE],
                                  check=False)

            # do not use return code, because this will not work when wrapped another time
            if EE.OE_HYBRID_HELP_IDENTIFICATION_STRING in result.stderr:
                return True
            return False
        except Exception as e:
            return False
