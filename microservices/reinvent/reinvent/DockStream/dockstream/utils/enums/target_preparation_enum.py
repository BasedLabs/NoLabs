
class TargetPreparationEnum:
    """This "Enum" serves to store all the strings used in parsing "target preparation" configurations."""

    TARGETPREP = "target_preparation"

    # header region
    # ---------
    HEADER = "header"
    ENVIRONMENT = "environment"
    ENVIRONMENT_EXPORT = "export"
    ENVIRONMENT_EXPORT_KEY = "key"
    ENVIRONMENT_EXPORT_VALUE = "value"

    INPUT_PATH = "input_path"

    # loggers
    # ---------
    LOGGING = "logging"
    LOGGING_ENABLED = "enabled"
    LOGGING_VERBOSITY = "verbosity"
    LOGGING_VERBOSITY_LOW = "low"
    LOGGING_VERBOSITY_HIGH = "high"
    LOGGING_LOGFILE = "logfile"

    # structure fixes
    # ---------
    FIX = "fixer"
    FIX_ENABLED = "enabled"
    FIX_STANDARDIZE = "standardize"
    FIX_MISSINGHEAVYATOMS = "fix_missing_heavy_atoms"
    FIX_MISSINGHYDROGENS = "fix_missing_hydrogens"
    FIX_MISSINGLOOPS = "fix_missing_loops"
    FIX_ADDWATERBOX = "add_water_box"
    FIX_REMOVEHETEROGENS = "remove_heterogens"
    FIX_PBDOUTPUTPATH = "fixed_pdb_path"

    # target preparation
    # ---------
    CAVITY = "cavity"
    CAVITY_METHOD = "method"

    CAVITY_METHOD_REFERENCE = "reference_ligand"
    CAVITY_REFERENCE_PATH = "reference_ligand_path"
    CAVITY_REFERENCE_FORMAT = "reference_ligand_format"
    CAVITY_REFERENCE_FORMAT_SDF = "SDF"
    CAVITY_REFERENCE_FORMAT_PDB = "PDB"

    # backend runs
    # ---------
    RUNS = "runs"

    RUNS_BACKEND = "backend"
    RUNS_BACKEND_RDOCK = "rDock"
    RUNS_BACKEND_OPENEYE = "OpenEye"
    RUNS_BACKEND_GOLD = "Gold"
    RUNS_BACKEND_AUTODOCKVINA = "AutoDockVina"

    RUNS_OUTPUT = "output"
    RUNS_PARAM = "parameters"
    RUNS_PARAM_PREFIX_EXECUTION = "prefix_execution"
    RUNS_PARAM_BINARY_LOCATION = "binary_location"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
