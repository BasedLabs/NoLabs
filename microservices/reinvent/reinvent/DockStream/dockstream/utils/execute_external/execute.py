import os
import abc
import subprocess
from shlex import quote


class ExecutorBase(metaclass=abc.ABCMeta):
    """Virtual base class for the general and program-specific executors."""

    def __init__(self, prefix_execution=None, binary_location=None):
        # if something needs to be attached to the execution string each time, store it here; if not, value is "None"
        self._prefix_execution = prefix_execution
        self._binary_location = binary_location

    @abc.abstractmethod
    def execute(self, command: str, arguments: list, check=True, location=None):
        # to avoid security issues, escape the arguments
        arguments = [quote(str(arg)) for arg in arguments]

        # check, if command (binary) is to be found at a specific location (rather than in $PATH)
        if self._binary_location is not None:
            command = os.path.join(self._binary_location, command)

        # check, if the something needs to be added before the execution of the "rDock" command
        if self._prefix_execution is not None:
            command = self._prefix_execution + " && " + command

        # execute; if "location" is set, change to this directory and execute there
        complete_command = [command + ' ' + ' '.join(str(e) for e in arguments)]
        old_cwd = os.getcwd()
        if location is not None:
            os.chdir(location)
        result = subprocess.run(complete_command,
                                check=check,                # force python to raise exception if anything goes wrong
                                universal_newlines=True,    # convert output to string (instead of byte array)
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                shell=True)
        os.chdir(old_cwd)
        return result

    @abc.abstractmethod
    def is_available(self):
        raise NotImplementedError("Overwrite this method in the child class.")


class Executor(ExecutorBase):
    """For execution of command-line programs that do not have any specific executor themselves."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(prefix_execution=prefix_execution, binary_location=binary_location)

    def execute(self, command: str, arguments: list, check=True, location=None):
        return super().execute(command=command,
                               arguments=arguments,
                               check=check,
                               location=location)

    def is_available(self):
        raise NotImplementedError("Cannot reliably check, whether a random program executes properly - do not use.")

