from dockstream.utils.enums.target_preparation_enum import TargetPreparationEnum
from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum
from dockstream.utils.enums.docking_enum import DockingConfigurationEnum, ResultKeywordsEnum


class GoldLigandPreparationEnum(LigandPreparationEnum):

    REMOVE_UNKNOWN_ATOMS = "remove_unknown_atoms"            # default: True; Whether or not to remove unknown atoms
    ASSIGN_BOND_TYPES = "assign_bond_types"                  # default: True; Whether or not to assign bond types
    STANDARDISE_BOND_TYPES = "standardise_bond_types"        # default: False; Whether or not to standardise bonds to CSD conventions
    ADD_HYDROGENS = "add_hydrogens"                          # default: True; Whether hydrogens need to be added
    PROTONATE = "protonate"                                  # default: True; Whether protonation rules need to be applied
    PROTONATION_RULES_FILE = "protonation_rules_file"        # default: None; Location of a file containing protonation rules

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class GoldTargetPreparationEnum(TargetPreparationEnum):

    OUTPUT_RECEPTORPATH = "receptor_path"

    CAVITY_REFERENCE_DISTANCE = "distance"

    CAVITY_METHOD_POINT = "point"
    CAVITY_POINT_ORIGIN = "origin"
    CAVITY_POINT_DISTANCE = "distance"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class GoldExecutablesEnum:
    """This "Enum" serves to store all the executables (and parameters) as strings available in the "Gold" module."""

    GOLD_AUTO = "gold_auto"
    GOLD_AUTO_HELP = "-h"
    GOLD_AUTO_HELP_IDENTIFICATION_STRING = "Usage: gold_auto"

    GOLD_AUTO_CONFIG_NAME = "gold_auto.config"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class GoldTargetKeywordEnum(GoldTargetPreparationEnum):

    VERSION = "version"
    CURRENT_VERSION = 1.0

    TARGET_PDB = "target_pdb"
    TARGET_PDB_FILENAME = "target_pdb_filename"
    REFERENCE_LIGAND = "reference_ligand"
    REFERENCE_LIGAND_FILENAME = "reference_ligand_filename"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class GoldDockingConfigurationEnum(DockingConfigurationEnum):

    RECEPTOR_PATHS = "receptor_paths"

    FITNESS_FUNCTION = "fitness_function"
    FITNESS_FUNCTION_GOLDSCORE = "goldscore"
    FITNESS_FUNCTION_CHEMSCORE = "chemscore"
    FITNESS_FUNCTION_ASP = "asp"
    FITNESS_FUNCTION_PLP = "plp"

    EARLY_TERMINATION = "early_termination"
    DIVERSE_SOLUTIONS = "diverse_solutions"

    NDOCKS = "ndocks"                             # number of docking attempts per ligand

    AUTOSCALE = "autoscale"                       # "very fast": 10
                                                  # "fast": 25
                                                  # "medium": 50
                                                  # "slow": 75
                                                  # "very slow": 100

    GOLD_RESPONSE_VALUE = "response_value"        # set to either "value" or "fitness" to use either as response
    GOLD_RESPONSE_VALUE_FITNESS = "fitness"
    GOLD_RESPONSE_VALUE_VALUE = "value"           # the "value" (and whether positive or negative values are "better")
                                                  # depend on the fitness function chosen

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class GoldResultKeywordsEnum(ResultKeywordsEnum):
    """This "Enum" serves to store all keywords for "Gold" result dictionaries."""

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class GoldOutputEnum:
    """This "Enum" serves to store all keywords that are used by the "Gold" module."""

    TAG = "tag"
    BEST = "best"
    DICT_FITNESS = {"chemscore": {TAG: "Gold.Chemscore.Fitness", BEST: "max"},
                    "plp": {TAG: "Gold.PLP.Fitness", BEST: "max"},
                    "goldscore": {TAG: "Gold.Goldscore.Fitness", BEST: "max"},
                    "asp": {TAG: "Gold.ASP.Fitness", BEST: "max"}}
    DICT_VALUE = {"asp": {TAG: "Gold.ASP.ASP", BEST: "max"},
                  "chemscore": {TAG: "Gold.Chemscore.DG", BEST: "min"},
                  "plp": {TAG: "Gold.PLP.PLP", BEST: "min"},
                  "goldscore": None}

    # The following values are "the-higher-the-better":
    # all fitness scores, Gold.ASP.ASP

    # The following values are "the-lower-the-better":
    # Gold.Chemscore.DG, Gold.PLP.PLP

    # For "goldscore", there is no "value" version.

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
