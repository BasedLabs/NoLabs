
class AnalysisConfigurationEnum:

    ANALYSIS = "analysis"

    # RMSD
    # ---------
    RMSD = "rmsd"
    RMSD_DATA = "data"
    RMSD_DATA_MOLECULES_PATH = "molecules_path"
    RMSD_DATA_NAME = "name"
    RMSD_OUTPUT = "output"
    RMSD_OUTPUT_SUMMARY_PATH = "summary_path"
    RMSD_OUTPUT_DETAILS_PATH = "details_path"
    RMSD_OUTPUT_HEATMAP_PATH = "heatmap_path"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class AnalysisInternalEnum:

    DATA_NAME = "name"
    MOLECULES = "molecules"
    FIRST_SET = "first_set"
    SECOND_SET = "second_set"
    LIST_RMSD_VALUES = "list_rmsd_values"
    MEAN_RMSD_VALUES = "mean_rmsd_values"
    SD_RMSD_VALUES = "sd_rmsd_values"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class AnalysisEnum:
    """this "Enum" defines all the strings required in the analysis script."""

    # Input Docking Data
    # --------------------------------------

    INPUT_DOCKING_DATA = "input_docking_data"
    DATA_PATH = "data_path"
    LIGAND_NUMBER = "ligand_number"
    DATA_METRIC = "data_metric"
    MAX_DATA_METRIC_BEST = "max_data_metric_best"
    DATA_THRESHOLDS = "data_thresholds"

    # Input Experimental Data
    # --------------------------------------

    INPUT_EXP_DATA = "input_exp_data"
    EXP_DATA_PATH = "exp_data_path"
    EXP_METRIC = "exp_metric"
    COMPARISON_SCORE = "comparison_score"
    MAX_EXP_METRIC_BEST = "max_exp_metric_best"
    EXP_THRESHOLDS = "exp_thresholds"

    # Input Binary Actives/Inactives Data
    # ---------------------------------------

    INPUT_ENRICHMENT_DATA = "input_enrichment_data"
    DATA_PATH_ACTIVES = "data_path_actives"
    DATA_PATH_INACTIVES = "data_path_inactives"
    ACTIVES_DATA_METRIC = "actives_data_metric"
    INACTIVES_DATA_METRIC = "inactives_data_metric"
    MAX_METRIC_BEST = "max_metric_best"
    ACTIVES = "Actives"
    INACTIVES = "Inactives"

    # Plots
    # ---------------------------------------

    PLOT_SETTINGS = "plot_settings"
    ENRICHMENT_ANALYSIS = "enrichment_analysis"
    PROC_OVERLAY = "pROC_overlay"

    # Output Folder
    # ---------------------------------------

    OUTPUT = "output"
    OUTPUT_PATH = "output_path"

    # Histogram Plot Parameters
    # ---------------------------------------

    HIST_TWO_COLOURS = ['green', 'blue']
    HIST_THREE_COLOURS = ['green', 'orange', 'blue']
    HIST_FOUR_COLOURS = ['green', 'orange', 'blue', 'cyan']

