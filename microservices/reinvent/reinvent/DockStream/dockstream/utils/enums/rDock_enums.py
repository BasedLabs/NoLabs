from dockstream.utils.enums.target_preparation_enum import TargetPreparationEnum
from dockstream.utils.enums.docking_enum import DockingConfigurationEnum, ResultKeywordsEnum


class rDockTargetPreparationEnum(TargetPreparationEnum):

    CAVITY_PRMFILE = "prm_file"
    CAVITY_METHOD_TWOSPHERES = "two_spheres"

    # strings that are replace in the PRM file
    STRING_RECEPTOR_MOL2_PATH = "<RECEPTOR_MOL2_FILE_ABSOLUTE_PATH>"
    STRING_REFERENCE_LIGAND_SDF_PATH = "<REFERENCE_SDF_FILE_ABSOLUTE_PATH>"

    RUNS_OUTPUT_DIRECTORY = "directory"

    PRM_DEFAULT_PATH = "dockstream/config/rDock/standard_reference_ligand.prm"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class rDockDockingConfigurationEnum(DockingConfigurationEnum):

    PARAMS_PRM_PATHS = "rbdock_prm_paths"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class rDockExecutablesEnum:
    """This "Enum" serves to store all the executables (and parameters) as strings available in the "rDock" module."""

    # environment variables
    # ---------
    RBT_ROOT = "RBT_ROOT"
    RBT_HOME = "RBT_HOME"

    # executable "rbdock" + parameters
    # ---------
    RBDOCK = "rbdock"
    RBDOCK_HELP = "--help"
    RBDOCK_HELP_IDENTIFICATION_STRING = "Usage: rbdock"
    RBDOCK_S = "-s"                           # random seed
    RBDOCK_S_DEFAULT = 42                     # default for random seed
    RBDOCK_T = "-T"                           # trace level for debugging
    RBDOCK_N = "-n"                           # number of runs
    RBDOCK_R = "-r"                           # receptor PRM file
    RBDOCK_O = "-o"                           # output file
    RBDOCK_P = "-p"                           # protocol file
    RBDOCK_P_DEFAULT = "dock.prm"             # default docking protocol, located in $RBT_ROOT/data/scripts
    RBDOCK_I = "-i"                           # input file

    # executable "rbcavity" + parameters
    # ---------
    RBCAVITY = "rbcavity"                     # generate a cavity from a receptor file
    RBCAVITY_R = "-r"                         # specifies path to input configuration PRM file
    RBCAVITY_D = "-d"                         # dump the grid in pymol
    RBCAVITY_WAS = "-was"                     # write the cavity out

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class rDockRbdockOutputEnum:
    """This "Enum" serves to store all keywords that are used by the "rbdock" executable."""

    # key values
    # ---------
    NAME = "Name"
    SCORE = "SCORE"

    # remaining values (for reference purposes only)
    # ---------
    CHROM_ZERO = "CHROM.0"
    CHROM_ONE = "CHROM.1"
    RI = "RI"

    RBT_CURRENT_DIRECTORY = "Rbt.Current_Directory"
    RBT_EXECUTABLE = "Rbt.Executable"
    RBT_LIBRARY = "Rbt.Library"
    RBT_PARAMETER_FILE = "Rbt.Parameter_File"
    RBT_RECEPTOR = "Rbt.Receptor"

    SCORE_INTER = "SCORE.INTER"
    SCORE_INTER_CONST = "SCORE.INTER.CONST"
    SCORE_INTER_POLAR = "SCORE.INTER.POLAR"
    SCORE_INTER_REPUL = "SCORE.INTER.REPUL"
    SCORE_INTER_ROT = "SCORE.INTER.ROT"
    SCORE_INTER_VDW = "SCORE.INTER.VDW"
    SCORE_INTER_NORM = "SCORE.INTER.norm"

    SCORE_INTRA = "SCORE.INTRA"
    SCORE_INTRA_DIHEDRAL = "SCORE.INTRA.DIHEDRAL"
    SCORE_INTRA_DIHEDRAL_ZERO = "SCORE.INTRA.DIHEDRAL.0"
    SCORE_INTRA_POLAR = "SCORE.INTRA.POLAR"
    SCORE_INTRA_POLAR_ZERO = "SCORE.INTRA.POLAR.0"
    SCORE_INTRA_REPUL = "SCORE.INTRA.REPUL"
    SCORE_INTRA_REPUL_ZERO = "SCORE.INTRA.REPUL.0"
    SCORE_INTRA_VDW = "SCORE.INTRA.VDW"
    SCORE_INTRA_VDW_ZERO = "SCORE.INTRA.VDW.0"
    SCORE_INTRA_NORM = "SCORE.INTRA.norm"

    SCORE_RESTR = "SCORE.RESTR"
    SCORE_RESTR_CAVITY = "SCORE.RESTR.CAVITY"
    SCORE_RESTR_NORM = "SCORE.RESTR.norm"

    SCORE_SYSTEM = "SCORE.SYSTEM"
    SCORE_SYSTEM_CONST = "SCORE.SYSTEM.CONST"
    SCORE_SYSTEM_DIHEDRAL = "SCORE.SYSTEM.DIHEDRAL"
    SCORE_SYSTEM_POLAR = "SCORE.SYSTEM.POLAR"
    SCORE_SYSTEM_REPUL = "SCORE.SYSTEM.REPUL"
    SCORE_SYSTEM_VDW = "SCORE.SYSTEM.VDW"
    SCORE_SYSTEM_NORM = "SCORE.SYSTEM.norm"

    SCORE_HEAVY = "SCORE.heavy"
    SCORE_NORM = "SCORE.norm"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class rDockRbcavityOutputEnum:
    """This "Enum" serves to store all keywords that are used by the "rbcavity" executable."""

    DOCKING_SITE = "DOCKING SITE"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class rDockResultKeywordsEnum(ResultKeywordsEnum):
    """This "Enum" serves to store all keywords for "rDock" result dictionaries."""

    # rDockTargetPreparator::specify_cavity() result dictionary
    # ---------
    SPECIFYCAVITY_BINARY_PATH = "binary_path"
    SPECIFYCAVITY_GRID_PATH = "grid_path"
    SPECIFYCAVITY_METADATA = "cavity_metadata"
    SPECIFYCAVITY_METADATA_TOTALVOLUME = "total_volume"
    SPECIFYCAVITY_METADATA_SIZEINPOINTS = "size_in_points"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
