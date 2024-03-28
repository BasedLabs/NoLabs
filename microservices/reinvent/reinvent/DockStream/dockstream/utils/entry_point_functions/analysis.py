import os
import sys
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from scipy.stats import spearmanr
from scipy.stats import kendalltau
import math
from dockstream.utils.enums.analysis_enums import AnalysisEnum

_AE = AnalysisEnum()


def correlation_analysis(parameters: dict, csv_file: str, analysis_results: dict):
    """This function takes a parameters dictionary provided by the user containing all required parameters and
       a results dictionary to be populated with statistical analysis metrics. A DataFrame is created containing
       the sorted (based on ligand number) and merged docking and experimental data. Correlation and ranked
       correlation metrics are calculated (Pearson, Spearman, Kendall) and stored in the results dictionary.
       The parameters dictionary (which may be modified if there are overlapping metric names from the user
       provided data), the results dictionary, and the DataFrame is returned for further plotting applications.

    :param parameters: a dictionary containing all the required parameters provided by the user
    :param csv_file: a single output file
    :param analysis_results: a dictionary to be populated with correlation and ranked correlation metrics
                           (Pearson, Spearman, Kendall
    :return: parameters dictionary, analysis_results dictionary, DataFrame containing the sorted (based on ligand
             number) and merged docking and experimental data
    """

    # generate the full path to the docking data
    if parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_PATH].endswith(".csv"):
        data_path = parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_PATH]

    else:
        data_path = os.path.join(parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_PATH], csv_file)
    # load the docking results
    docking_data = pd.read_csv(data_path)
    # load the experimental comparison data
    experimental_data = pd.read_csv(parameters[_AE.INPUT_EXP_DATA][_AE.EXP_DATA_PATH])

    # drop all columns except ligand_number and scores in the docking data
    docking_data.drop(docking_data.columns.difference([_AE.LIGAND_NUMBER, parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_METRIC]]), 1, inplace=True)
    # drop all columns except ligand_number and scores in the experimental data
    # if the user specified docking and exp data have the same metric name (ex. "score"), rename
    # the exp_data_metric to comparison_score to ensure the DataFrame is unambiguous
    if parameters[_AE.INPUT_EXP_DATA][_AE.EXP_METRIC] == parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_METRIC]:
        parameters[_AE.INPUT_EXP_DATA][_AE.EXP_METRIC] = _AE.COMPARISON_SCORE
        experimental_data.rename(columns={parameters[_AE.INPUT_EXP_DATA][_AE.EXP_METRIC]: _AE.COMPARISON_SCORE}, inplace=True)
    # the below code block is a safeguard for when the user wants to run analysis for multiple csv files
    # where both the docking and exp data share the same metric name (ex. "score")
    if parameters[_AE.INPUT_EXP_DATA][_AE.EXP_METRIC] == _AE.COMPARISON_SCORE:
        experimental_data.rename(columns={parameters[_AE.INPUT_EXP_DATA][_AE.EXP_METRIC]: _AE.COMPARISON_SCORE}, inplace=True)

    experimental_data.drop(experimental_data.columns.difference([_AE.LIGAND_NUMBER, parameters[_AE.INPUT_EXP_DATA][_AE.EXP_METRIC]]), 1, inplace=True)
    # merge the docking results and experimental binding data based on their ligand number identifiers.
    # This is to ensure that data pertaining to the correct ligands are compared irrespective of the
    # order in which the user provides the experimental data
    comparison_df = docking_data.merge(experimental_data, on=_AE.LIGAND_NUMBER)

    # extract the docking scores and experimental potency parameter columns
    # also correct for sign discrepancies (ex. if higher docking scores are better or lower docking scores are better)
    docking_scores, experimental_data = corr_sign_corrector(parameters, comparison_df)
    # fit a linear model using least squares
    model = LinearRegression().fit(docking_scores, experimental_data)
    coeff_determination = model.score(docking_scores, experimental_data)
    # perform Spearman ranked correlation analysis
    spearman_correlation = spearmanr(docking_scores, experimental_data)
    # perform Kendall ranked correlation analysis
    kendall_correlation = kendalltau(docking_scores, experimental_data)
    # store the Coefficient of Determination, Spearman, and Kendall quantities in a dictionary
    analysis_results[csv_file] = {"coeff_determination": coeff_determination, "Spearman_coeff":
                                  spearman_correlation[0], "Kendall_coeff": kendall_correlation[0]}

    return parameters, analysis_results, comparison_df


def corr_sign_corrector(parameters: dict, comparison_df: pd.DataFrame) -> np.array:
    """this function takes a parameters dictionary provided by the user containing all required parameters and DataFrame
       containing the sorted (based on ligand number) and merged docking and experimental data. The docking data and
       experimental data are extracted and transformed based on the user specified max_data_metric_best and
       max_exp_metric_best booleans. This ensures the computed statistics (e.g. Spearman Coefficient) follows the
       conventional signs for easier interpretability (i.e. "higher" docking scores correlating to "higher" exp. metric
       yields a positive Spearman Coefficient)

    :param parameters: a dictionary containing all the required parameters provided by the user
    :param comparison_df: a DataFrame containing the sorted (based on ligand number) and merged docking
                          and experimental data
    :return: 2 numpy arrays containing the extracted and modified (if applicable) docking and experimental data
    """

    # interpret the user max_data_metric_best boolean
    max_data_metric_best = parameters[_AE.INPUT_DOCKING_DATA][_AE.MAX_DATA_METRIC_BEST]
    if max_data_metric_best.lower() == "true":
        max_data_metric_best = True
    elif max_data_metric_best.lower() == "false":
        max_data_metric_best = False
    else:
        sys.exit("'max_data_metric_best' value error: 'enrichment_analysis' is set to 'False'. Assuming either "
                 "Correlation Analysis or Thresholds Analysis is desired, please set 'max_data_metric_best' to"
                 "either 'True' or 'False'. Exiting script.")

    # interpret the user max_exp_metric_best boolean
    max_exp_metric_best = parameters[_AE.INPUT_EXP_DATA][_AE.MAX_EXP_METRIC_BEST]
    if max_exp_metric_best.lower() == "true":
        max_exp_metric_best = True
    elif max_exp_metric_best.lower() == "false":
        max_exp_metric_best = False
    else:
        sys.exit("'max_exp_metric_best' value error: 'enrichment_analysis' is set to 'False'. Assuming either "
                 "Correlation Analysis or Thresholds Analysis is desired, please set 'max_exp_metric_best' to"
                 "either 'True' or 'False'. Exiting script.")

    if not max_data_metric_best and not max_exp_metric_best:
        docking_scores = comparison_df[parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_METRIC]].apply(lambda x: x * -1).to_numpy().reshape(-1, 1)
        experimental_data = comparison_df[parameters[_AE.INPUT_EXP_DATA][_AE.EXP_METRIC]].apply(lambda x: x * -1).to_numpy()

    elif max_data_metric_best and not max_exp_metric_best:
        docking_scores = comparison_df[parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_METRIC]].to_numpy().reshape(-1, 1)
        experimental_data = comparison_df[parameters[_AE.INPUT_EXP_DATA][_AE.EXP_METRIC]].apply(lambda x: x * -1).to_numpy()

    elif not max_data_metric_best and max_exp_metric_best:
        docking_scores = comparison_df[parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_METRIC]].apply(lambda x: x * -1).to_numpy().reshape(-1, 1)
        experimental_data = comparison_df[parameters[_AE.INPUT_EXP_DATA][_AE.EXP_METRIC]].to_numpy()

    elif max_data_metric_best and max_exp_metric_best:
        docking_scores = comparison_df[parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_METRIC]].to_numpy().reshape(-1, 1)
        experimental_data = comparison_df[parameters[_AE.INPUT_EXP_DATA][_AE.EXP_METRIC]].to_numpy()

    return docking_scores, experimental_data


def binary_data_classification(parameters: dict, comparison_df: pd.DataFrame, dock_thresh: float, exp_thresh: float):
    """this function takes a DataFrame containing the sorted (based on ligand number) and merged docking and
       experimental data and the parameters dictionary provided by the user containing all required parameters.
       The docking and experimental data points are compared to the user specified docking/exp thresholds to create
       a binary classification encompassing "active" and "inactive". The binarized docking and experimental data is
       returned for further plotting applications

    :param parameters: a dictionary containing all the required parameters provided by the user
    :param comparison_df: a DataFrame containing the sorted (based on ligand number) and merged docking
                          and experimental data
    :param dock_thresh: a float representing the user provided docking threshold
    :param exp_thresh: a float representing the user provided experimental threshold
    :return: 2 lists containing the binarized docking and experimental data
    """

    # interpret the user max_data_metric_best boolean
    max_data_metric_best = parameters[_AE.INPUT_DOCKING_DATA][_AE.MAX_DATA_METRIC_BEST]
    if max_data_metric_best.lower() == "true":
        max_data_metric_best = True
    elif max_data_metric_best.lower() == "false":
        max_data_metric_best = False
    else:
        sys.exit("'max_data_metric_best' value error: set to either 'True' or 'False'. Exiting the script.")

    # interpret the user max_exp_metric_best boolean
    max_exp_metric_best = parameters[_AE.INPUT_EXP_DATA][_AE.MAX_EXP_METRIC_BEST]
    if max_exp_metric_best.lower() == "true":
        max_exp_metric_best = True
    elif max_exp_metric_best.lower() == "false":
        max_exp_metric_best = False
    else:
        sys.exit("'max_exp_metric_best' value error: set to either 'True' or 'False'. Exiting the script.")

    # load the docking scores and experimental data into lists to allow easy classification into pairs
    # encompassing true positive, false positives, true negatives, and false negatives
    docking_scores = comparison_df[parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_METRIC]].tolist()
    exp_data = comparison_df[parameters[_AE.INPUT_EXP_DATA][_AE.EXP_METRIC]].tolist()

    # initialize lists to be populated with binarized docking and experimental data
    binary_docking_scores = []
    binary_exp_data = []

    if not max_data_metric_best and not max_exp_metric_best:
        for val in docking_scores:
            if val <= dock_thresh:
                binary_docking_scores.append(1)
            else:
                binary_docking_scores.append(0)

        for val in exp_data:
            if val <= exp_thresh:
                binary_exp_data.append(1)
            else:
                binary_exp_data.append(0)

    elif max_data_metric_best and not max_exp_metric_best:
        for val in docking_scores:
            if val >= dock_thresh:
                binary_docking_scores.append(1)
            else:
                binary_docking_scores.append(0)

        for val in exp_data:
            if val <= exp_thresh:
                binary_exp_data.append(1)
            else:
                binary_exp_data.append(0)

    elif not max_data_metric_best and max_exp_metric_best:
        for val in docking_scores:
            if val <= dock_thresh:
                binary_docking_scores.append(1)
            else:
                binary_docking_scores.append(0)

        for val in exp_data:
            if val >= exp_thresh:
                binary_exp_data.append(1)
            else:
                binary_exp_data.append(0)

    elif max_data_metric_best and max_exp_metric_best:
        for val in docking_scores:
            if val >= dock_thresh:
                binary_docking_scores.append(1)
            else:
                binary_docking_scores.append(0)

        for val in exp_data:
            if val >= exp_thresh:
                binary_exp_data.append(1)
            else:
                binary_exp_data.append(0)

    return binary_docking_scores, binary_exp_data


def data_classification(parameters: dict, comparison_df: pd.DataFrame, dock_thresh: float, exp_thresh: float):
    """this function takes a DataFrame containing the sorted (based on ligand number) and merged docking and
       experimental data and the parameters dictionary provided by the user containing all required parameters.
       The docking and experimental data points are compared to the user specified docking/exp thresholds to
       create a quaternary classification encompassing true_pos, false_pos, true_neg, false_neg. The classified
       data points are returned for further plotting applications

    :param parameters: a dictionary containing all the required parameters provided by the user
    :param comparison_df: a DataFrame containing the sorted (based on ligand number) and merged docking
                          and experimental data
    :return: 4 lists containing the classified true_pos, false_pos, true_neg, false_neg data points
    """

    # interpret the user max_data_metric_best boolean
    max_data_metric_best = parameters[_AE.INPUT_DOCKING_DATA][_AE.MAX_DATA_METRIC_BEST]
    if max_data_metric_best.lower() == "true":
        max_data_metric_best = True
    elif max_data_metric_best.lower() == "false":
        max_data_metric_best = False
    else:
        sys.exit("'max_data_metric_best' value error: set to either 'True' or 'False'. Exiting the script.")

    # interpret the user max_exp_metric_best boolean
    max_exp_metric_best = parameters[_AE.INPUT_EXP_DATA][_AE.MAX_EXP_METRIC_BEST]
    if max_exp_metric_best.lower() == "true":
        max_exp_metric_best = True
    elif max_exp_metric_best.lower() == "false":
        max_exp_metric_best = False
    else:
        sys.exit("'max_exp_metric_best' value error: set to either 'True' or 'False'. Exiting the script.")

    # load the docking scores and experimental data into lists to allow easy classification into pairs
    # encompassing true positive, false positives, true negatives, and false negatives
    docking_scores = comparison_df[parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_METRIC]].tolist()
    exp_data = comparison_df[parameters[_AE.INPUT_EXP_DATA][_AE.EXP_METRIC]].tolist()

    if not max_data_metric_best and not max_exp_metric_best:
        true_pos = [[dock, exp] for dock, exp in zip(docking_scores, exp_data) if
                    dock <= dock_thresh and exp <= exp_thresh]
        false_pos = [[dock, exp] for dock, exp in zip(docking_scores, exp_data) if
                     dock <= dock_thresh and exp > exp_thresh]
        true_neg = [[dock, exp] for dock, exp in zip(docking_scores, exp_data) if
                    dock > dock_thresh and exp > exp_thresh]
        false_neg = [[dock, exp] for dock, exp in zip(docking_scores, exp_data) if
                     dock > dock_thresh and exp <= exp_thresh]

    elif max_data_metric_best and not max_exp_metric_best:
        true_pos = [[dock, exp] for dock, exp in zip(docking_scores, exp_data) if
                    dock >= dock_thresh and exp <= exp_thresh]
        false_pos = [[dock, exp] for dock, exp in zip(docking_scores, exp_data) if
                     dock >= dock_thresh and exp > exp_thresh]
        true_neg = [[dock, exp] for dock, exp in zip(docking_scores, exp_data) if
                    dock < dock_thresh and exp > exp_thresh]
        false_neg = [[dock, exp] for dock, exp in zip(docking_scores, exp_data) if
                     dock < dock_thresh and exp <= exp_thresh]

    elif not max_data_metric_best and max_exp_metric_best:
        true_pos = [[dock, exp] for dock, exp in zip(docking_scores, exp_data) if
                    dock <= dock_thresh and exp >= exp_thresh]
        false_pos = [[dock, exp] for dock, exp in zip(docking_scores, exp_data) if
                     dock <= dock_thresh and exp < exp_thresh]
        true_neg = [[dock, exp] for dock, exp in zip(docking_scores, exp_data) if
                    dock > dock_thresh and exp < exp_thresh]
        false_neg = [[dock, exp] for dock, exp in zip(docking_scores, exp_data) if
                     dock > dock_thresh and exp >= exp_thresh]

    elif max_data_metric_best and max_exp_metric_best:
        true_pos = [[dock, exp] for dock, exp in zip(docking_scores, exp_data) if
                    dock >= dock_thresh and exp >= exp_thresh]
        false_pos = [[dock, exp] for dock, exp in zip(docking_scores, exp_data) if
                     dock >= dock_thresh and exp < exp_thresh]
        true_neg = [[dock, exp] for dock, exp in zip(docking_scores, exp_data) if
                    dock < dock_thresh and exp < exp_thresh]
        false_neg = [[dock, exp] for dock, exp in zip(docking_scores, exp_data) if
                     dock < dock_thresh and exp >= exp_thresh]

    return true_pos, false_pos, true_neg, false_neg


def plot_settings(true_pos: list, false_pos: list, true_neg: list, false_neg: list, file_name: str,
                  dock_thresh: float, exp_thresh: float):
    """this function takes the classified points (true_pos, false_pos, true_neg, false_neg) and determines if any of
    the lists are empty. Based on which data type(s) (ex. true_neg), if applicable, is(are) empty, plot settings
    including the number of labels required and the number of colours required are determined and returned for
    histogram plotting.

    :param true_pos: a list containing the true_pos data
    :param false_pos: a list containing the false_pos data
    :param true_neg: a list containing the true_neg data
    :param false_neg: a list containing the false_neg data
    :param file_name: the name of the current DockStream output file being analyzed,
                      used to print out logging messages for the user
    :param dock_thresh: a float representing the user provided docking threshold
    :param exp_thresh: a float representing the user provided experimental threshold
    :return: plot labels and plot colours
    """

    # initialize dictionaries to keep track of which data sets are empty
    data_name_tracker = {0: 'true pos', 1: 'false pos', 2: 'true neg', 3: 'false neg'}
    # list to keep track of which data sets (names) are populated
    populated_labels = []
    # index counter to iterate through dictionary
    tracker_index = 0

    for classification in (true_pos, false_pos, true_neg, false_neg):
        if len(classification) != 0:
            populated_labels.append(data_name_tracker[tracker_index])
            tracker_index += 1
        else:
            print(f"Note: During histogram data processing for {file_name}, the specified docking threshold: "
                  f"{dock_thresh} and exp threshold: {exp_thresh} failed to yield any {data_name_tracker[tracker_index]} docked/exp pairs")
            tracker_index += 1

    # based on the number of classified data sets that are not empty, set marker and colours to preset settings
    if len(populated_labels) == 2:
        hist_colours = _AE.HIST_TWO_COLOURS
    if len(populated_labels) == 3:
        hist_colours = _AE.HIST_THREE_COLOURS
    if len(populated_labels) == 4:
        hist_colours = _AE.HIST_FOUR_COLOURS

    return populated_labels, hist_colours


def pROC_curve_datapoints(parameters: dict, actives_data: pd.DataFrame, inactives_data: pd.DataFrame) -> list:
    """this function takes the parameters dictionary provided by the user containing all required parameters.
       The actives and inactives data are iterated to generate (x, y) data points for pROC curve construction
       and pROC AUC calculation. The data points are returned for plotting applications

    :param parameters: a dictionary containing all the required parameters provided by the user
    :param actives_data: a DataFrame containing the actives data
    :param inactives_data: a DataFrame containing the inactives data
    :return: 4 lists containing the false positive rates (pROC curve x variable), true positive rates (pROC curve y variable),
             pROC curve (x, y) data points for the random classifier, and all the false positive rates required for pROC AUC
    """
    # interpret the enrichment analysis max_metric_best boolean
    max_metric_best = parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.MAX_METRIC_BEST]
    if max_metric_best.lower() == "true":
        max_metric_best = True
    elif max_metric_best.lower() == "false":
        max_metric_best = False
    else:
        sys.exit("'max value best' value error: set to either 'True' or 'False'. Exiting the script without generating pROC curves.")

    if max_metric_best:
        actives_data = [(score, "active") for score in actives_data[parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.ACTIVES_DATA_METRIC]]]
        inactives_data = [(score, "inactive") for score in inactives_data[parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.INACTIVES_DATA_METRIC]]]
    elif not max_metric_best:
        actives_data = [(-score, "active") for score in actives_data[parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.ACTIVES_DATA_METRIC]]]
        inactives_data = [(-score, "inactive") for score in inactives_data[parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.INACTIVES_DATA_METRIC]]]

    all_data = actives_data + inactives_data
    all_data.sort(reverse=True)

    # initialize a lower bound as the log function is undefined at 0
    lower_bound = 1 / len(inactives_data)
    TPR, FPR, rand_selection, FPR_AUC = [lower_bound], [lower_bound], [lower_bound], []
    num_actives, num_inactives = 0, 0

    for score, tag in all_data:
        if tag == "active":
            num_actives += 1
        if tag == "inactive":
            num_inactives += 1

        TPR_value = num_actives / len(actives_data)
        FPR_value = num_inactives / len(inactives_data)

        # cap the lower bound value
        if FPR_value < lower_bound:
            FPR_value = lower_bound

        TPR.append(TPR_value)
        FPR.append(FPR_value)
        # FPR = TPR for a random classifier so just append FPR values and use it for both x and y axes during plotting
        rand_selection.append(FPR_value)

        # pROC AUC calculation requires the FPR values corresponding to each time an active is recovered
        if tag == "active":
            FPR_AUC.append(FPR_value)

    return TPR, FPR, rand_selection, FPR_AUC


def pROC_AUC(FPR_AUC: list) -> float:
    """this function takes all the false positive rates each corresponding to when an active ligand was recovered.
       This list is generated by pROC_curve_datapoints. The pROC AUC is calculated and returned

    :param FPR_AUC: a list containing all the false positive rates corresponding to when an active ligand was recovered
    :return: pROC AUC value
    """

    pROC_AUC_sum = 0
    for value in FPR_AUC:
        pROC_AUC_sum += math.log10((1 / value))

    pROC_AUC = round((pROC_AUC_sum / len(FPR_AUC)), 3)

    return pROC_AUC


def enrichment_factor(parameters: dict, actives_data: pd.DataFrame, inactives_data: pd.DataFrame, EF_values: dict, csv_name="Docking Experiment EF 5%") -> dict:
    """this function calculates the enrichment factor (EF) for the top 5% of ligands based on docking score.
       A dictionary containing the stored EF 5% value is returned

    :param parameters: a dictionary containing all the required parameters provided by the user
    :param actives_data: a DataFrame containing the actives data
    :param inactives_data: a DataFrame containing the inactives data
    :param EF_values: a dictionary to store all calculated EF 5% values
    :param csv_name: the name of the current data file being analyzed. This parameter is only relevant during
                     batch enrichment analysis
    :return: a dictionary containing the stored EF 5% value(s)
    """

    # interpret the enrichment analysis max_metric_best boolean
    max_metric_best = parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.MAX_METRIC_BEST]
    if max_metric_best.lower() == "true":
        max_metric_best = True
    elif max_metric_best.lower() == "false":
        max_metric_best = False
    else:
        sys.exit("'max_metric_best' value error: set to either 'True' or 'False'. Exiting the script.")

    if max_metric_best:
        actives_data = [(score, "active") for score in actives_data[parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.ACTIVES_DATA_METRIC]]]
        inactives_data = [(score, "inactive") for score in inactives_data[parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.INACTIVES_DATA_METRIC]]]
    elif not max_metric_best:
        actives_data = [(-score, "active") for score in actives_data[parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.ACTIVES_DATA_METRIC]]]
        inactives_data = [(-score, "inactive") for score in inactives_data[parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.INACTIVES_DATA_METRIC]]]

    all_data = actives_data + inactives_data
    all_data.sort(reverse=True)

    num_actives_ef = 0
    # count how many active ligands are in the top 5% of sorted ligands
    for score, tag in all_data[:int(len(all_data) * 0.05)]:
        if tag == "active":
            num_actives_ef += 1

    file_name = csv_name.replace(".csv", "")
    EF_values[file_name] = round((num_actives_ef / (len(actives_data) * 0.05)), 2)

    return EF_values











