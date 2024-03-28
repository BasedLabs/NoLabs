from dockstream.utils.enums.target_preparation_enum import TargetPreparationEnum
from dockstream.utils.enums.docking_enum import DockingConfigurationEnum

from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum


class LigprepLigandPreparationEnum(LigandPreparationEnum):

    USE_EPIK = "use_epik"
    TARGET_PH = "target_pH"
    PH_TOLERANCE = "pH_tolerance"
    CHIRALITY = "chirality"
    MAX_NUMBER_STEREOISOMERS = "max_number_stereoisomers"
    FORCE_FIELD = "force_field"
    FORCE_FIELD_OPLS_2005 = "OPLS_2005"
    FORCE_FIELD_OPLS3e = "OPLS3e"
    FILTER_FILE = "filter_file"
    COMMAND_LINE_PARAMETERS = "command_line_parameters"

    # parallelization
    # ---------
    PARALLELIZATION = "parallelization"
    PARALLELIZATION_NUMBER_CORES = "number_cores"
    PARALLELIZATION_MAXCOMPOUNDSPERSUBJOB = "max_compounds_per_subjob"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class SchrodingerDockingConfigurationEnum(DockingConfigurationEnum):

    GLIDE_FLAGS = "glide_flags"
    GLIDE_KEYWORDS = "glide_keywords"

    GLIDE_INPUTBLOCK_COMMASEPARATED = ["CONSTRAINT_GROUP"]     # define list of block keys which are to have commas
    GLIDE_INPUTBLOCK_VALUEQUOTED = ["FEATURE"]                 # define list of block keys, where values are to be put into double quotation marks

    GLIDE_TIME_LIMIT_PER_COMPOUND = "time_limit_per_compound"

    GLIDE_TG = "token_guard"
    GLIDE_TG_PREFIX_EXECUTION = "prefix_execution"
    GLIDE_TG_BINARY_LOCATION = "binary_location"
    GLIDE_TG_TOKEN_POOLS = "token_pools"
    GLIDE_TG_WAIT_INTERVAL = "wait_interval_seconds"
    GLIDE_TG_WAIT_LIMIT = "wait_limit_seconds"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class SchrodingerExecutablesEnum:
    """This "Enum" serves to store all the executables (and parameters and special strings) as strings available in the "Schrodinger" module."""

    # executable "sdconvert" + parameters
    # ---------
    SDCONVERT = "sdconvert"
    SDCONVERT_CALL = "$SCHRODINGER/utilities/sdconvert"
    SDCONVERT_HELP = ""
    SDCONVERT_HELP_IDENTIFICATION_STRING = "mae : Maestro format"
    SDCONVERT_A = "-a"                                               # append structures to the output file
    SDCONVERT_I = "-i"                                               # input; note that the format is directly appended (e.g. "-isd")
    SDCONVERT_O = "-o"                                               # output; note that the format is directly appended (e.g. "-omae")
    SDCONVERT_FORMAT_SD = "sd"                                       # MDL SDfile format
    SDCONVERT_FORMAT_MM = "mm"                                       # MacroModel (.dat) format
    SDCONVERT_FORMAT_MAE = "mae"                                     # Maestro format
    SDCONVERT_TITLE = "-title"                                       # define SD property <prop> as the source of the Maestro title
    SDCONVERT_NOSTEREO = "-nostereo"                                 # do not record the atom parity info from the input file
    SDCONVERT_NOAROM = "-noarom"                                     # do not convert aromatic type 4 bonds to single and double bonds (which is the Maestro convention)

    # executable "ligprep" + parameters
    # ---------
    LIGPREP = "ligprep"
    LIGPREP_CALL = "$SCHRODINGER/ligprep"
    LIGPREP_INPUT_ISMI = "-ismi"                                     # SMI input followed by <path> (alternatives: "-icsv", "-imae" and "-isd")
    LIGPREP_OUTPUT_OSD = "-osd"                                      # SD(F) output followed by <path> (alternative: "-omae")
    LIGPREP_INP_CONFIG = "-inp"                                      # not used in DockStream, but would be an option to feed parameters from configuration file
    LIGPREP_EPIK = "-epik"                                           # Use "Epik" for ionization and tautomerization (Recommended)
    LIGPREP_PH = "-ph"                                               # Effective / target pH; followed by <number> (use 7.0 as default)
    LIGPREP_PHT = "-pht"                                             # pH tolerance for generated structures; followed by <number> (use 2.0 as default)
    LIGPREP_AC = "-ac"                                               # Do not respect existing chirality properties and do not respect chiralities from the input geometry. Generate stereoisomers for all chiral centers up to
                                                                     # the number permitted (specified using the -s option). This is equivalent to "Generate all combinations" in the Ligand Preparation user interface. Default
                                                                     # behavior is to respect only explicitly indicated chiralities.
    LIGPREP_F = "-f"                                                 # Filter structures via LigFilter using specifications from the file provided. Default: do not filter.
    LIGPREP_G = "-g"                                                 # Respect chiralities from input geometry when generating stereoisomers.
    LIGPREP_S = "-s"                                                 # Generate up to this <number> stereoisomers per input structure. (Default: 32).
    LIGPREP_BFF = "-bff"                                             # Force-field to be used for the final geometry optimization (either 14 or 16, which refers to OPLS_2005 and
                                                                     # OPLS3e respectively. Default: 14
    LIGPREP_FF_OPLS_2005 = "14"                                      # Default force-field
    LIGPREP_FF_OPLS3e = "16"                                         # Alternative force-field
    LIGPREP_NJOBS = "-NJOBS"                                         # Divide the overall job into NJOBS subjobs. Set to 1 by default.
    LIGPREP_NSTRUCTS = "-NSTRUCTS"                                   # Divide the overall job into subjobs with no more than NSTRUCTS structures. Set to 1 by default.
    LIGPREP_HOST = "-HOST"                                           # Run the job on <hostname> remotely on the indicated host entry.
    LIGPREP_HOST_LOCALHOST = "localhost"                             # Default value for the run.
    LIGPREP_WAIT = "-WAIT"                                           # Do not return a prompt until the job completes.

    # executable "licadmin" + parameters
    # ---------
    LICADMIN = "licadmin"
    LICADMIN_STAT = "STAT"                                           # returns the list of tokens used / available

    # executable "glide" + parameters
    # note, that you can get the full list of parameters with "$SCHRODINGER/glide -k"
    # ---------
    GLIDE = "glide"
    GLIDE_CALL = "$SCHRODINGER/glide"
    GLIDE_HELP = "-h"
    GLIDE_HELP_IDENTIFICATION_STRING = "positional arguments:"
    GLIDE_WAIT = "-WAIT"
    GLIDE_OVERWRITE = "-OVERWRITE"                                   # Remove previous job files before running.
    GLIDE_NJOBS = "-NJOBS"                                           # Divide the overall job into NJOBS subjobs.
    GLIDE_HOST = "-HOST"                                             # Run job remotely on the indicated host entry.
    GLIDE_TMPLAUNCHDIR = "-TMPLAUNCHDIR"                             # WARNING: does not seem to be supported (any longer?) - probably "-NOLOCAL" now?
    GLIDE_ATTACHED = "-ATTACHED"                                     # WARNING: does not seem to be supported (any longer?)
    GLIDE_AMIDE_MODE = "AMIDE_MODE"                                  # amide bond rotation behavior: "fixed", "free", "penal", "trans", "gen[eralized]"
    GLIDE_EXPANDED_SAMPLING = "EXPANDED_SAMPLING"                    # bypass elimination of poses in rough scoring stage (useful for fragment docking)
    GLIDE_GRIDFILE = "GRIDFILE"                                      # path to grid (.grd or .zip) file
    GLIDE_LIGANDFILE = "LIGANDFILE"                                  # Glide docking ligands file name
    GLIDE_NENHANCED_SAMPLING = "NENHANCED_SAMPLING"                  # expand size of the Glide funnel by N times to process poses from N confgen runs with minor perturbations to the input ligand coordinates
    GLIDE_POSE_OUTTYPE = "POSE_OUTTYPE"                              # format for file containing docked poses: "poseviewer" for _pv.mae output; "ligandlib" for _lib.mae; similarly "poseviewer_sd" and "ligandlib_sd" for sdf output; "phase_subset" for bypassing _lib or _pv in favor of a Phase subset file.
    GLIDE_POSE_OUTTYPE_LIGANDLIB = "ligandlib_sd"                    # the only supported type for now; is enforce in the backend
    GLIDE_POSES_PER_LIG = "POSES_PER_LIG"                            # maximum number of poses to report per each input ligand
    GLIDE_POSTDOCK_NPOSE = "POSTDOCK_NPOSE"                          # maximum number of best-by-Emodel poses to submit to post-docking minimization
    GLIDE_POSTDOCKSTRAIN = "POSTDOCKSTRAIN"                          # include strain correction in post-docking score
    GLIDE_PRECISION = "PRECISION"                                    # glide docking precision ("SP", "XP" or "HTVS")
    GLIDE_REWARD_INTRA_HBONDS = "REWARD_INTRA_HBONDS"                # reward formation of intramolecular hydrogen bonds in the ligand
    GLIDE_USE_CONS = "USE_CONS"
    GLIDE_NREQUIRED_CONS = "NREQUIRED_CONS"
    GLIDE_LOG_SUCCESS_STRING = "glide_sort command succeeded"        # if any of these string is present in the logfile associated with a subjob, all went well
    GLIDE_LOG_FINISHED_STRINGS = {"Exiting Glide"}
    GLIDE_LOG_FAIL_STRINGS = {"*** Error in",                        # if any of these strings is present in the logfile associated with a subjob, there was an issue resulting in the complete failure of the execution
                              "Glide cannot recover from this signal and will now abort.",
                              "======= Backtrace: ========="}
                              #"Glide: FATAL mmlewis error"}

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class SchrodingerOutputEnum:
    """This "Enum" serves to store all keywords that are used by the "Schrodinger" module."""

    # Glide
    # ---------
    GLIDE_DOCKING_SCORE = "r_i_docking_score"                    # the docking score (including "Epik" corrections")
    GLIDE_GSCORE = "r_i_glide_gscore"                            # the "docking score" without "Epik" corrections
    GLIDE_SOURCE_FILE_INDEX = "i_m_source_file_index"            # the index of the ligand in the input file (starting with '1')

    GLIDE_SDF_DEFAULT_EXTENSION = "_lib.sdfgz"
    GLIDE_LOG = ".log"
    GLIDE_SDF = ".sdf"

    # Ligprep
    # ---------
    LIGPREP_LOG = ".log"
    LIGPREP_VARIANTS = "s_lp_Variant"                            # the SDF tag with <identifier>-# (where # is the number of the enumeration starting with '1')
    LIGPREP_TAUTOMER_PROBABILITY = "r_lp_tautomer_probability"   # number from 0 to 1 (sums up to 1 over all variants)

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class SchrodingerTargetPreparationEnum(TargetPreparationEnum):


    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
