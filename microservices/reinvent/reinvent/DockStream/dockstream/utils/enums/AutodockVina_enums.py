from dockstream.utils.enums.target_preparation_enum import TargetPreparationEnum
from dockstream.utils.enums.docking_enum import DockingConfigurationEnum, ResultKeywordsEnum


class AutodockVinaDockingConfigurationEnum(DockingConfigurationEnum):

    ADV_RECEPTOR_PTBQT_PATH = "receptor_pdbqt_path"
    ADV_SEED = "seed"
    ADV_SEARCH_SPACE = "search_space"
    ADV_SEARCH_SPACE_CENTER_X = "--center_x"
    ADV_SEARCH_SPACE_CENTER_Y = "--center_y"
    ADV_SEARCH_SPACE_CENTER_Z = "--center_z"
    ADV_SEARCH_SPACE_SIZE_X = "--size_x"
    ADV_SEARCH_SPACE_SIZE_Y = "--size_y"
    ADV_SEARCH_SPACE_SIZE_Z = "--size_z"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class AutodockVinaExecutablesEnum:

    # executable "vina" + parameters
    # ---------
    VINA = "vina"
    VINA_CALL = "vina"                                 # the binary call
    VINA_HELP = "--help"                               # display usage summary
    VINA_HELP_ADVANCED = "--help_advanced"             # display usage summary (with all options)
    VINA_VERSION = "--version"                         # diplay program version
    VINA_VERSION_IDENTIFICATION_STRING = "AutoDock Vina v1.2.5"   # string, which needs to be present in help output in
                                                                 # order to assume "AutoDock Vina" can be properly used
    VINA_CONFIGURATION = "--config"                    # path to configuration file, where options below can be put

    # input
    VINA_RECEPTOR = "--receptor"                       # rigid part of the receptor (PDBQT)
    VINA_LIGAND = "--ligand"                           # ligand (PDBQT); only one at a time
    VINA_FLEX = "--flex"                               # flexible side chains, if any (PDBQT)

    # search space
    VINA_CENTER_X = "--center_x"                       # X coordinate of the center
    VINA_CENTER_Y = "--center_y"                       # Y coordinate of the center
    VINA_CENTER_Z = "--center_z"                       # Z coordinate of the center
    VINA_SIZE_X = "--size_x"                           # size in the X dimension (Angstroms)
    VINA_SIZE_Y = "--size_y"                           # size in the X dimension (Angstroms)
    VINA_SIZE_Z = "--size_z"                           # size in the X dimension (Angstroms)

    # output
    VINA_OUT = "--out"                                 # output models (PDBQT), the default is chosen based on the
                                                       # ligand file name
    VINA_LOG = "--log"                                 # optionally, write log file

    # advanced options
    VINA_SCORE_ONLY = "--score_only"                   # score only - search space can be omitted
    VINA_LOCAL_ONLY = "--local_only"                   # do local search only
    VINA_RANDOMIZE_ONLY = "--randomize_only"           # randomize input, attempting to avoid clashes
    VINA_WEIGHT_GAUSS1 = "--weight_gauss1"             # gauss_1 weight (default: -0.035579)
    VINA_WEIGHT_GAUSS2 = "--weight_gauss2"             # gauss_2 weight (default: -0.005156)
    VINA_WEIGHT_REPULSION = "--weight_repulsion"       # repulsion weight (default: 0.84024500000000002)
    VINA_WEIGHT_HYDROPHOBIC = "--weight_hydrophobic"   # hydrophobic weight (-0.035069000000000003)
    VINA_WEIGHT_HYDROGEN = "--weight_hydrogen"         # hydrogen bond weight (-0.58743900000000004)
    VINA_WEIGHT_ROT = "--weight_rot"                   # N_rot weight (default: 0.058459999999999998)

    # miscellaneous (optional)
    VINA_CPU = "--cpu"                                 # the number of CPUs to use (the default is to try to detect
                                                       # the number of CPUs or, failing that, use 1)
    VINA_SEED = "--seed"                               # explicit random seed
    VINA_EXHAUSTIVENESS = "--exhaustiveness"           # exhaustiveness of the global search (roughly proportional
                                                       # to time): 1+ (default: 8)
    VINA_NUM_MODES = "--num_modes"                     # maximum number of binding modes to generate (default: 9)
    VINA_ENERGY_RANGE = "--energy_range"               # maximum energy difference between the best binding mode and the
                                                       # worst one displayed [kcal/mol] (default: 3)

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class AutodockVinaOutputEnum:

    ADV_PDBQT = ".pdbqt"

    # the score is part of a tag in the PDBQT -> SDF translated output (tag "REMARK"), which looks like that:
    # < REMARK >
    # VINA RESULT: -9.1 0.000 0.000
    # Name = /tmp/tmpjssiy8z4.pdb
    # ...

    # Note, that the three values are: affinity [kcal/mol] | dist from best mode (rmsd l.b.) | rmsd (u. b.)
    REMARK_TAG = "REMARK"
    RESULT_LINE_IDENTIFIER = "VINA RESULT"
    RESULT_LINE_POS_SCORE = 2
    RESULT_LINE_POS_RMSDTOBEST_LB = 3
    RESULT_LINE_POS_RMSDTOBEST_UB = 4

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class AutodockTargetPreparationEnum(TargetPreparationEnum):

    ADV_PDBQT = ".pdbqt"
    RECEPTOR_PATH = "receptor_path"
    PH = "pH"
    EXTRACT_BOX = "extract_box"
    EXTRACT_BOX_REFERENCE_LIGAND_PATH = "reference_ligand_path"
    EXTRACT_BOX_REFERENCE_LIGAND_FORMAT = "reference_ligand_format"
    EXTRACT_BOX_REFERENCE_LIGAND_FORMAT_PDB = "PDB"
    EXTRACT_BOX_REFERENCE_LIGAND_FORMAT_SDF = "SDF"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class AutodockResultKeywordsEnum(ResultKeywordsEnum):
    """This "Enum" serves to store all keywords for "AutoDock Vina" result strings."""

    SDF_TAG_SCORE = "SCORE"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
