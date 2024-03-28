
class LigandPreparationEnum:
    """This "Enum" serves to store all the strings used in parsing "DockStream" configurations concerning ligands."""

    LIGAND_PREPARATION = "ligand_preparation"

    # the embedding (3D coordinate generation)
    EMBEDDING_POOLS = "embedding_pools"
    POOLID = "pool_id"
    MOLECULES = "molecules"
    PARAMS = "parameters"

    # input parameters
    INPUT = "input"
    INPUT_STANDARDIZE_SMILES = "standardize_smiles"
    INPUT_PATH = "input_path"
    INPUT_TYPE = "type"
    INPUT_TYPE_CONSOLE = "CONSOLE"
    INPUT_TYPE_LIST = "LIST"
    INPUT_TYPE_SMI = "SMI"
    INPUT_TYPE_CSV = "CSV"
    INPUT_TYPE_SDF = "SDF"
    PREFIX_EXECUTION = "prefix_execution"
    BINARY_LOCATION = "binary_location"
    INITIALIZATION_MODE = "initialization_mode"
    INITIALIZATION_MODE_ORDER = "order"
    INITIALIZATION_MODE_AZDOCK = "dockstream"

    # CSV input
    INPUT_CSV_DELIMITER = "delimiter"
    INPUT_CSV_DELIMITER_DEFAULT = ','
    INPUT_CSV_COLUMNS = "columns"
    INPUT_CSV_COLNAME_SMILES = "smiles"
    INPUT_CSV_COLNAME_NAMES = "names"

    # SDF input
    INPUT_SDF_TAGS = "tags"
    INPUT_SDF_TAGNAME_NAMES = "names"

    # output parameters
    OUTPUT = "output"
    OUTPUT_CONFORMERPATH = "conformer_path"
    OUTPUT_FORMAT = "format"
    OUTPUT_FORMAT_SDF = "SDF"
    OUTPUT_FORMAT_MOL2 = "MOL2"

    # TautEnum can be used to prepare the smiles (input)
    # ---------
    USE_TAUT_ENUM = "use_taut_enum"
    TAUT_ENUM_PREFIX_EXECUTION = "prefix_execution"
    TAUT_ENUM_ENUMERATE_PROTONATION = "enumerate_protonation"
    TAUT_ENUM_BINARY_LOCATION = "binary_location"

    # the different types of embedding
    TYPE = "type"
    TYPE_RDKIT = "RDkit"
    TYPE_OPENEYE = "OpenEye"
    TYPE_CORINA = "Corina"
    TYPE_LIGPREP = "Ligprep"
    TYPE_GOLD = "Gold"
    TYPE_OMEGA = "Omega"

    # structural alignment to reference ("internal" method)
    # ---------
    ALIGN = "align"
    ALIGN_MODE = "mode"
    ALIGN_MODE_INTERNAL = "internal"
    ALIGN_REFERENCE_PATHS = "reference_paths"
    ALIGN_REFERENCE_FORMAT = "reference_format"
    ALIGN_REFERENCE_FORMAT_SDF = "SDF"
    ALIGN_REFERENCE_FORMAT_PDB = "PDB"
    ALIGN_MINIMUM_SUBSTRUCTURE_RATIO = "minimum_substructure_ratio"
    ALIGN_COMPLETE_RINGS_ONLY = "complete_rings_only"

    # what to do in cases, where no alignment can be made, define what to do
    ALIGN_FAIL_ACTION = "fail_action"
    ALIGN_FAIL_DISCARD = "discard"
    ALIGN_FAIL_KEEP = "keep"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
