from dockstream.utils.enums.target_preparation_enum import TargetPreparationEnum
from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum
from dockstream.utils.enums.docking_enum import DockingConfigurationEnum, ResultKeywordsEnum


class OpenEyeLigandPreparationEnum(LigandPreparationEnum):

    # align using OpenEye's template version, which is set at the receptor building stage
    ALIGN_MODE_OPENEYERECEPTOR = "OpenEye_receptor"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class OpenEyeTargetPreparationEnum(TargetPreparationEnum):

    OUTPUT_RECEPTORPATH = "receptor_path"

    CAVITY_METHOD_BOX = "box"
    CAVITY_BOX_LIMITS = "limits"
    CAVITY_METHOD_HINT = "hint"
    CAVITY_HINT_COORDINATES = "coordinates"


    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class OpenEyeDockingConfigurationEnum(DockingConfigurationEnum):

    RECEPTOR_PATHS = "receptor_paths"

    # scoring functions
    # ---------
    SCORING = "scoring"
    SCORING_INVALID_VALUE = 16777215
    # McGann2003: shape-based scoring function that favours poses that complement the active site well, ignoring any
    # chemical interactions; good choice to ensure shape-complementarity
    SCORING_SHAPEGAUSS = "Shapegauss"
    # Verkhivker2002: Piecewise Linear Potential uses both shape and hydrogen bond complementarity; in the implementation
    # used, it also includes metal-based interactions
    SCORING_PLP = "PLP"
    # Eldridge1997: includes lipophilic, H-bonds, metals, clashes, rotatable bonds
    SCORING_CHEMSCORE = "Chemscore"
    # the Chemgauss-scoring functions use Gaussian smoothed potentials to measure complementarity; includes shape,
    # H-bonds between ligand and protein, H-bonds with implicit solvent and metal interactions; version 4 is an
    # improvement in terms of H-bonding
    SCORING_CHEMGAUSS3 = "Chemgauss3"
    SCORING_CHEMGAUSS4 = "Chemgauss4"
    SCORING_HYBRID1 = "Hybrid1"
    SCORING_HYBRID2 = "Hybrid2"

    # resolution (specifies search resolution during exhaustive search and local optimization as well as the number
    # of poses passed from the exhaustive step to the optimization step
    # ---------
    RESOLUTION = "resolution"
    RESOLUTION_INVALID_VALUE = 16777215
    RESOLUTION_HIGH = "High"              # 1000 poses passed
    RESOLUTION_STANDARD = "Standard"      # 100 poses passed
    RESOLUTION_LOW = "Low"                # 100 poses passed

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class OpenEyeResultKeywordsEnum(ResultKeywordsEnum):
    """This "Enum" serves to store all keywords for "OpenEye" result dictionaries."""



    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
