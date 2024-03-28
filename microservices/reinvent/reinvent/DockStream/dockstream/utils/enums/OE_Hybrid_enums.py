from dockstream.utils.enums.target_preparation_enum import TargetPreparationEnum
from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum
from dockstream.utils.enums.docking_enum import DockingConfigurationEnum, ResultKeywordsEnum


class OpenEyeHybridLigandPreparationEnum(LigandPreparationEnum):

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


class OpenEyeHybridTargetPreparationEnum(TargetPreparationEnum):

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


class OpenEyeHybridExecutablesEnum:
    """This "Enum" serves to store all the executables (and parameters) as strings available in the "OpenEye Hybrid" module"""

    OE_HYBRID_MODULE_LOAD = "module load oedocking"

    HYBRID = "hybrid"
    OE_HYBRID_HELP_SIMPLE = "--help simple"   # returns simple list of parameters
    OE_HYBRID_HELP_ALL = "--help all"    # returns complete list of parameters
    OE_HYBRID_HELP_DEFAULTS = "--help defaults"   # returns the default values for all parameters
    OE_HYBRID_HELP_HTML = "--help html"   # creates an html help file for OE Hybrid
    OE_HYBRID_HELP_VERSIONS = "--help versions"   # lists toolkits and versions used for OE Hybrid
    OE_HYBRID_HELP_IDENTIFICATION_STRING = "To cite HYBRID"   # string to identify whether OE Hybrid is available for a docking job
    # required parameters
    # -------------------
    RECEPTOR = "-receptor"   # required: receptor file for docking. Must contain bound ligand
    DBASE = "-dbase"   # required: database; ligands to dock. Ligand must be in 3D format
    # optional parameters
    # -------------------
    PARAM = "-param"   # text file containing parameters for docking job
    MOLNAMES = "-molnames"   # text file containing molecule names corresponding to ligands in the dbase parameter. Only ligands with matched names will be docked
    DOCK_RESOLUTION = "-dock_resolution"   # controls docking resolution, default: standard
    DOCKED_MOLECULE_FILE = "-docked_molecule_file"   # file to write the docked molecules to, default: docked.oeb.gz
    UNDOCKED_MOLECULES_FILE = "-undocked_molecule_file"   # file to write the unsuccessfully docked molecules to, default: undocked.oeb.gz
    SCORE_FILE = "-score_file"   # file to write the docked molecule names and scores to, default: score.txt
    REPORT_FILE = "-report_file"   # file to write the text report of the docking run to, default: report.txt
    SETTINGS_FILE = "-settings_file"   # file to write the settings used for the docking run, default: settings.param
    STATUS_FILE = "-status_file"   # file to write the status of the docking run which is updated and overwritten every few seconds, default: status.txt
    HITLIST_SIZE = "-hitlist_size"   # parameter controls the number of top scoring molecules to be outputted. Excess will be discarded.
                                     # "0" denotes "serial mode" where all molecules will be outputted unsorted, default: 500
    NUM_POSES = "-num_poses"   # parameter specifies number of docked poses to output for each molecules
    SCORE_TAG = "-score_tag"   # parameter specifies tag to use when storing molecule scores, default: HYBRID Chemgauss4 Score
    ANNOTATE_SCORES = "-annotate_scores"   # parameter specifies whether to add VIDA (OpenEye's molecular visualization program) score annotations to processed molecules, defautlt: false
    SAVE_COMPONENT_SCORES = "-save_component_scores"   # parameter specifies whether individual components of the total score for each pose is saved to the score file
    NO_EXTRA_OUTPUT_FILES = "-no_extra_output_files"   # parameter controls output files from docking run. if "true", the output is only the docked molecule file, default: false
    NO_DOTS = "-no_dots"   # parameter specifies whether a dot/"x" is written for standard error/failed docking molecules, default: false
    PREFIX = "-prefix"   # parameter specifies prefix to use for all output files, default: hybrid

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class OpenEyeHybridDockingConfigurationEnum(DockingConfigurationEnum):

    RECEPTOR_PATHS = "receptor_paths"
    NUMBER_POSES = "number_poses"

    # scoring functions
    # ---------
    SCORING = "scoring"
    SCORING_INVALID_VALUE = 16777215
    # McGann2003: shape-based scoring function that favours poses that complement the active site well, ignoring any
    # chemical interactions; good choice to ensure shape-complementarity
    SCORING_SHAPEGAUSS = "Shapegauss"
    # Verkhivker2000: Piecewise Linear Potential uses both shape and hydrogen bond complementarity; in the implementation
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


class OpenEyeHybridResultKeywordsEnum(ResultKeywordsEnum):
    """This "Enum" serves to store all keywords for "OpenEye Hybrid" result dictionaries"""

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class OpenEyeHybridOutputKeywordsEnum:
    """This "Enum" serves to store all "OpenEye Hybrid" output enums"""

    # there are 6 default output files by default and in some cases (e.g. docked poses), the file
    # extension provided denotes the output format

    SCORE = "HYBRID Chemgauss4 score"

    DOCKED_MOLECULES_SDF_OUTPUT = "docked_molecules.sdf"   # specifying sdf here removes need for OpenBabel conversion of oeb to sdf
    UNDOCKED_MOLECULES_SDF_OUTPUT = "undocked_molecules.sdf"   # specifying sdf here removes need for OpenBabel conversion of oeb to sdf
    SCORE_FILE_OUTPUT = "score.txt"
    REPORT_FILE_OUTPUT = "report.txt"
    SETTINGS_FILE_OUTPUT = "settings.param"
    STATUS_FILE_OUTPUT = "status.txt"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
