
class DockingConfigurationEnum:
    """This "Enum" serves to store all the strings used in parsing "DockStream" configurations."""

    DOCKING = "docking"

    # header region
    # ---------
    HEADER = "header"
    ENVIRONMENT = "environment"
    ENVIRONMENT_EXPORT = "export"
    ENVIRONMENT_EXPORT_KEY = "key"
    ENVIRONMENT_EXPORT_VALUE = "value"
    LOGGING = "logging"
    LOGGING_VERBOSITY = "verbosity"
    LOGGING_LOGFILE = "logfile"

    DOCKING_RUNS = "docking_runs"

    RUN_ID = "run_id"
    INPUT_POOLS = "input_pools"
    PARAMS = "parameters"
    PARAMS_PREFIX_EXECUTION = "prefix_execution"
    PARAMS_BINARY_LOCATION = "binary_location"

    # parallelization
    # ---------
    PARALLELIZATION = "parallelization"
    PARALLELIZATION_NUMBER_CORES = "number_cores"
    PARALLELIZATION_MAXCOMPOUNDSPERSUBJOB = "max_compounds_per_subjob"

    # the different backend types
    # ---------
    BACKEND = "backend"
    BACKEND_RDOCK = "rDock"
    BACKEND_OPENEYE = "OpenEye"
    BACKEND_OPENEYEHYBRID = "Hybrid"
    BACKEND_AUTODOCKVINA = "AutoDockVina"
    BACKEND_GOLD = "Gold"
    BACKEND_GLIDE = "Glide"

    # structural alignment to reference
    # ---------
    ALIGN = "align"
    ALIGN_REFERENCE_PATHS = "reference_paths"
    ALIGN_MINIMUM_SUBSTRUCTURE_RATIO = "minimum_substructure_ratio"
    ALIGN_COMPLETE_RINGS_ONLY = "complete_rings_only"

    # what to do in cases, where no alignment can be made, define what to do
    ALIGN_FAIL_ACTION = "fail_action"
    ALIGN_FAIL_DISCARD = "discard"
    ALIGN_FAIL_KEEP = "keep"

    # output
    OUTPUT = "output"
    OUTPUT_POSES = "poses"
    OUTPUT_POSES_OVERWRITE = "overwrite"
    OUTPUT_POSES_PATH = "poses_path"
    OUTPUT_SCORES = "scores"
    OUTPUT_SCORES_OVERWRITE = "overwrite"
    OUTPUT_SCORES_PATH = "scores_path"
    OUTPUT_MODE = "mode"
    OUTPUT_MODE_ALL = "all"
    OUTPUT_MODE_BESTPERLIGAND = "best_per_ligand"
    OUTPUT_MODE_BESTPERENUMERATION = "best_per_enumeration"

    # number poses returned from docking (top X poses)
    NUMBER_POSES = "number_poses"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class ResultKeywordsEnum:
    """This "Enum" serves to store all keywords for result dictionaries."""

    # ResultParser::get_result() result dataframe
    # ---------
    DF_LIGAND_NAME = "name"
    DF_LIGAND_NAME_MOLECULE = ""
    DF_LIGAND_NAME_CONFORMER = ""
    DF_LIGAND_NUMBER = "ligand_number"
    DF_LIGAND_ENUMERATION = "enumeration"
    DF_CONFORMER = "conformer_number"
    DF_SCORE = "score"
    DF_SMILES = "smiles"
    DF_LOWEST_CONFORMER = "lowest_conformer"
    AGGREGATE_BEST = "best"
    AGGREGATE_AVERAGE = "average"

    # fixed values
    # ---------
    FIXED_VALUE_NA = "NA"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
