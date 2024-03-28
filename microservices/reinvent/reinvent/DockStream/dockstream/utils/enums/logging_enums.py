
class LoggingConfigEnum:
    """This "Enum" serves to store all paths and keywords used to configure the loggers."""

    # set levels (for now, they match to the "logging" default ones)
    DEBUG = "debug"
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    EXCEPTION = "exception"

    # paths to the configuration JSONs that are shipped with DockStream
    PATH_CONFIG_DEFAULT = "dockstream/config/logging/default.json"
    PATH_CONFIG_VERBOSE = "dockstream/config/logging/verbose.json"
    PATH_CONFIG_DEBUG = "dockstream/config/logging/debug.json"

    # high-level loggers defined in the configurations
    LOGGER_INTERFACE = "command_line_interface"
    LOGGER_TARGET_PREPARATION = "target_preparation"
    LOGGER_LIGAND_PREPARATION = "ligand_preparation"
    LOGGER_DOCKING = "docking"
    LOGGER_BLANK = "blank"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
