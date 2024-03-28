import os
import sys
import errno
import json
import argparse
import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix
from sklearn.metrics import matthews_corrcoef
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sb

from dockstream.utils.enums.analysis_enums import AnalysisEnum
from dockstream.utils.entry_point_functions import analysis

_AE = AnalysisEnum()


def scat_plot(parameters: dict, csv_file: str, comparison_df: pd.DataFrame):
    """this function takes a DataFrame containing the sorted (based on ligand number) and merged docking and
       experimental data to generate a scatter plot. The parameters dictionary provided by the user containing
       all required parameters

    :param parameters: a dictionary containing all the required parameteres provided by the user
    :param csv_file: a single DockStream output file
    :param comparison_df: a DataFrame containing the sorted (based on ligand number) and merged docking
                          and experimental data
    """

    # load the docking scores and experimental data into arrays
    docking_scores = comparison_df[parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_METRIC]].to_numpy().reshape(-1, 1)
    experimental_data = comparison_df[parameters[_AE.INPUT_EXP_DATA][_AE.EXP_METRIC]].to_numpy()

    file_name = csv_file.replace(".csv", "")
    sb.set()
    plt.scatter(docking_scores, experimental_data, color="g")
    plt.title(file_name + " Scatter Plot")
    plt.xlabel(parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_METRIC])
    plt.ylabel(parameters[_AE.INPUT_EXP_DATA][_AE.EXP_METRIC])
    output_path = os.path.join(parameters[_AE.OUTPUT][_AE.OUTPUT_PATH], file_name)
    plt.savefig(output_path + "_scatplot.png", format="png", dpi=300)
    plt.figure()


def binary_matrix(parameters: dict, csv_file: str, comparison_df: pd.DataFrame) -> dict:
    """this function takes a DataFrame containing the sorted (based on ligand number) and merged docking and
       experimental data and the parameters dictionary provided by the user containing all required parameters.
       The docking and experimental data points are compared to the user specified docking/exp thresholds to create
       a binary classification encompassing "active" and "inactive". The binarized data is used to generate a confusion
       matrix and a dictionary containing the Matthews Correlation Coefficients (MCCs) for all trials is returned

    :param parameters: a dictionary containing all the required parameters provided by the user
    :param csv_file: a single DockStream output file
    :param comparison_df: a DataFrame containing the sorted (based on ligand number) and merged docking
                          and experimental data
    :return: dictionary containing the calculated Matthew Correlation Coefficient (MCC)
    """

    # initialize a dictionary to store calculated Matthew Correlation Coefficients (MCCs)
    matt_coeffs = {}

    # loop through all the user provided docking/exp threshold pairs (they are allowed to provide >1)
    # and perform confusion matrix classification analysis
    for threshold in zip(parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_THRESHOLDS], parameters[_AE.INPUT_EXP_DATA][_AE.EXP_THRESHOLDS]):
        dock_thresh, exp_thresh = threshold[0], threshold[1]

        # call a helper function to binarize the docking and experimental data
        binary_docking_scores, binary_exp_data = analysis.binary_data_classification(parameters, comparison_df, dock_thresh, exp_thresh)

        # calculate and store the Matthews Correlation Coefficient (MCC)
        matt_coeffs[f"MCC_for_{csv_file}_{str(threshold)}"] = matthews_corrcoef(binary_exp_data, binary_docking_scores)

        true_pos, false_pos, true_neg, false_neg = 0, 0, 0, 0

        # count the number of true_pos, false_pos, true_neg, and false_neg points
        for exp_point, data_point in zip(binary_exp_data, binary_docking_scores):
            if exp_point == 1 and data_point == 1:
                true_pos += 1
            if exp_point == 0 and data_point == 1:
                false_pos += 1
            if exp_point == 0 and data_point == 0:
                true_neg += 1
            if exp_point == 1 and data_point == 0:
                false_neg += 1

        total_data_points = true_pos + false_pos + true_neg + false_neg

        labels = np.asarray([[f"{true_neg} ({round(true_neg / total_data_points * 100, 2)}%) \n Dock Threshold = {dock_thresh} \n Exp. Threshold = {exp_thresh}",
                              f"{false_pos} ({round(false_pos / total_data_points * 100, 2)}%)"],
                             [f"{false_neg} ({round(false_neg / total_data_points * 100, 2)}%)",
                              f"{true_pos} ({round(true_pos / total_data_points * 100, 2)}%)"]])

        # generate the confusion matrix
        file_name = csv_file.replace(".csv", "")
        matrix = confusion_matrix(binary_exp_data, binary_docking_scores)
        sb.heatmap(matrix, annot=labels, fmt='', cmap='BuGn', cbar=False)
        plt.title(file_name + " Confusion Matrix")
        plt.xlabel("Binarized " + parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_METRIC])
        plt.ylabel("Binarized " + parameters[_AE.INPUT_EXP_DATA][_AE.EXP_METRIC])
        output_path = os.path.join(parameters[_AE.OUTPUT][_AE.OUTPUT_PATH], file_name)
        plt.savefig(output_path + f"_matrix_{str(threshold)}.png", format="png", dpi=300)
        plt.figure()

    # return the dictionary storing the Matthew Correlation Coefficients (MCCs) for all user provided threshold pairs
    return matt_coeffs


def histogram(parameters: dict, csv_file: str, comparison_df: pd.DataFrame):
    """this function takes a DataFrame containing the sorted (based on ligand number) and merged docking and
       experimental data and the parameters dictionary provided by the user containing all required parameters.
       The docking and experimental data points are compared to the user specified docking/exp thresholds to create
       a quaternary classification of true positives, false positives, true negatives, false negatives. Histograms
       are constructed to visualize the data distribution/investigate general separation of the data.

    :param parameters: a dictionary containing all the required parameters provided by the user
    :param csv_file: a single DockStream output file
    :param comparison_df: a DataFrame containing the sorted (based on ligand number) and merged docking
                          and experimental data
    """

    # generate the output path
    file_name = csv_file.replace(".csv", "")
    output_path = os.path.join(parameters[_AE.OUTPUT][_AE.OUTPUT_PATH], file_name)

    # loop through all the user provided docking/exp threshold pairs (they are allowed to provide >1)
    # and perform classification
    for threshold in zip(parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_THRESHOLDS], parameters[_AE.INPUT_EXP_DATA][_AE.EXP_THRESHOLDS]):
        dock_thresh, exp_thresh = threshold[0], threshold[1]

        # call helper function to perform data classification
        true_pos, false_pos, true_neg, false_neg = analysis.data_classification(parameters, comparison_df, dock_thresh, exp_thresh)

        # call helper function to generate plot settings for histogram applications
        populated_labels, hist_colours = analysis.plot_settings(true_pos, false_pos, true_neg, false_neg, file_name, dock_thresh, exp_thresh)

        # extract the docking scores corresponding to the classified points
        true_pos_points = [pair[0] for pair in true_pos]
        false_pos_points = [pair[0] for pair in false_pos]
        true_neg_points = [pair[0] for pair in true_neg]
        false_neg_points = [pair[0] for pair in false_neg]

        # remove any classifications that are not populated (ex. no true_neg points given the provided docking/exp thresholds)
        hist_points = []
        for points in (true_pos_points, false_pos_points, true_neg_points, false_neg_points):
            if len(points) != 0:
                hist_points.append(points)

        plt.hist(hist_points, color=hist_colours, label=populated_labels, bins=20, stacked=True, alpha=0.7)

        plt.title(file_name + " Histogram")
        plt.xlabel(parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_METRIC])
        plt.ylabel("Absolute Counts")
        plt.legend()
        plt.savefig(output_path + f"_histogram_{str(threshold)}.png", format="png", dpi=300)
        plt.figure()

        # if either true_pos or true_neg is empty (no classified points), then a pooled histogram is not applicable
        if len(true_pos) == 0 or len(true_neg) == 0:
            if len(true_pos) == 0 and len(true_neg) == 0:
                print(f"Note: During histogram data processing for {file_name}, the specified docking threshold: "
                      f"{dock_thresh} and exp threshold: {exp_thresh} failed to yield any True Pos and True Neg docked/exp pairs.\nPooled histograms cannot be generated.")
                continue
            elif len(true_pos) == 0:
                print(f"Note: During histogram data processing for {file_name}, the specified docking threshold: "
                      f"{dock_thresh} and exp threshold: {exp_thresh} failed to yield any True Pos docked/exp pairs.\nPooled histograms cannot be generated.")
                continue
            elif len(true_neg) == 0:
                print(
                    f"Note: During histogram data processing for {file_name}, the specified docking threshold: "
                    f"{dock_thresh} and exp threshold: {exp_thresh} failed to yield any True Neg docked/exp pairs.\nPooled histograms cannot be generated.")
                continue

        pooled_pos_points = true_pos_points + false_neg_points
        pooled_neg_points = true_neg_points + false_pos_points

        plt.hist([pooled_pos_points, pooled_neg_points], color=["g", "b"],
                 label=["Pooled Pos", "Pooled Neg"], bins=20, stacked=True, alpha=0.7)

        plt.title(file_name + " Pooled Density Plot")
        plt.xlabel(parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_METRIC])
        plt.ylabel("Absolute Counts")
        plt.legend(["Pooled Pos", "Pooled Neg"])
        plt.savefig(output_path + f"_pooled_histogram_{str(threshold)}.png", format="png", dpi=300)
        plt.figure()


def enrichment_analysis(parameters: dict, pROC_AUC_values: dict, EF_values: dict):
    """this function is only called if the user specifies for enrichment analysis The parameters dictionary provided
       by the user containing all required parameters. The pROC_AUC_values and EF_values dictionaries stores all the
       calculated pROC AUC and EF 5% values, respectively. Helper functions are called to generate enrichment
       histograms, box plots, and pROC curves

    :param parameters: a dictionary containing all the required parameters provided by the user
    :param pROC_AUC_values: a dictionary to store all calculated pROC AUC values
    :param EF_values: a dictionary to store all calculated EF 5% values
    """

    # check data_path_actives is valid:
    if not os.path.isdir(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_ACTIVES]) and not os.path.isfile(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_ACTIVES]):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_ACTIVES])
    # check data_path_inactives is valid:
    if not os.path.isdir(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_INACTIVES]) and not os.path.isfile(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_INACTIVES]):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_INACTIVES])

    # single DockStream run enrichment analysis
    if os.path.isfile(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_ACTIVES]) and os.path.isfile(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_INACTIVES]):
        # load the actives and inactives data
        actives_data = pd.read_csv(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_ACTIVES])
        inactives_data = pd.read_csv(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_INACTIVES])

        actives_df = pd.DataFrame({_AE.ACTIVES: actives_data[parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.ACTIVES_DATA_METRIC]].tolist()})
        inactives_df = pd.DataFrame({_AE.INACTIVES: inactives_data[parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.INACTIVES_DATA_METRIC]].tolist()})

        # generate enrichment histogram
        enrichment_histogram(parameters, actives_df, inactives_df)

        # generate enrichment boxplot
        enrichment_boxplot(parameters, actives_df, inactives_df)

        # generate pROC curve data points
        TPR, FPR, rand_selection, FPR_AUC = analysis.pROC_curve_datapoints(parameters, actives_data, inactives_data)

        # generate pROC curve and calculated pROC AUC
        pROC_AUC_values = pROC_curve(parameters, TPR, FPR, rand_selection, FPR_AUC, pROC_AUC_values)

        # calculate EF 5%
        EF_values = analysis.enrichment_factor(parameters, actives_data, inactives_data, EF_values)

        # save the pROC AUC values in json format for readability
        json_output_path = os.path.join(parameters[_AE.OUTPUT][_AE.OUTPUT_PATH], "pROC_AUC_values.json")
        with open(json_output_path, "w+") as f:
            json.dump(pROC_AUC_values, f, indent=2)

        # save the EF 5% values in json format for readability
        json_output_path = os.path.join(parameters[_AE.OUTPUT][_AE.OUTPUT_PATH], "EF_5%_values.json")
        with open(json_output_path, "w+") as f:
            json.dump(EF_values, f, indent=2)

        sys.exit("Enrichment histogram, boxplot, and pROC curve constructed. Exiting script now.")

    # batch DockStream runs enrichment analysis
    elif os.path.isdir(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_ACTIVES]) and os.path.isdir(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_INACTIVES]):
        # if user specifies to generate only an overlay pROC curve. The script exits if the below function call executes
        pROC_overlay_bool = parameters[_AE.PLOT_SETTINGS][_AE.PROC_OVERLAY]
        if pROC_overlay_bool.lower() == "true":
            pROC_overlay_bool = True
        elif pROC_overlay_bool.lower() == "false":
            pROC_overlay_bool = False
        else:
            sys.exit("'pROC_overlay' value error: set to either 'True' or 'False'. Exiting the script.")

        if pROC_overlay_bool:
            pROC_overlay_curve(parameters, pROC_AUC_values, EF_values)
        # proceed with normal Enrichment Analysis
        for csv_name in [csv for csv in os.listdir(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_ACTIVES]) if csv.endswith(".csv")]:
            # try (and bypass any error) to generate the enrichment plots
            # errors can occur in cases where the same experiment file name is not in both the actives and inactives folder
            try:
                # load the actives and inactives data
                actives_data = pd.read_csv(os.path.join(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_ACTIVES], csv_name))
                inactives_data = pd.read_csv(os.path.join(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_INACTIVES], csv_name))

                actives_df = pd.DataFrame({_AE.ACTIVES: actives_data[parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.ACTIVES_DATA_METRIC]].tolist()})
                inactives_df = pd.DataFrame({_AE.INACTIVES: inactives_data[parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.INACTIVES_DATA_METRIC]].tolist()})

                # generate enrichment histogram
                enrichment_histogram(parameters, actives_df, inactives_df, csv_name)

                # generate enrichment boxplot
                enrichment_boxplot(parameters, actives_df, inactives_df, csv_name)

                # generate pROC curve data points
                TPR, FPR, rand_selection, FPR_AUC = analysis.pROC_curve_datapoints(parameters, actives_data, inactives_data)

                # generate pROC curve
                pROC_curve(parameters, TPR, FPR, rand_selection, FPR_AUC, pROC_AUC_values, csv_name)

                # calculate EF 5%
                EF_values = analysis.enrichment_factor(parameters, actives_data, inactives_data, EF_values, csv_name)

            except:
                continue

    # save the pROC AUC values in json format for readability
    json_output_path = os.path.join(parameters[_AE.OUTPUT][_AE.OUTPUT_PATH], "pROC_AUC_values.json")
    with open(json_output_path, "w+") as f:
        json.dump(pROC_AUC_values, f, indent=2)

    # save the EF 5% values in json format for readability
    json_output_path = os.path.join(parameters[_AE.OUTPUT][_AE.OUTPUT_PATH], "EF_5%_values.json")
    with open(json_output_path, "w+") as f:
        json.dump(EF_values, f, indent=2)

    sys.exit("Batch enrichment histograms, boxplots, and pROC curves constructed. Exiting script now.")


def enrichment_histogram(parameters: dict, actives_df: pd.DataFrame, inactives_df: pd.DataFrame, csv_name=None):
    """this function takes the parameters dictionary provided by the user containing all required parameters.
       The actives and inactives data is used to construct an enrichment histogram

    :param parameters: a dictionary containing all the required parameters provided by the user
    :param actives_df: a DataFrame containing the actives data
    :param inactives_df: a DataFrame containing the inactives data
    :param csv_name: the name of the current data file being analyzed. This parameter is only relevant during
                     batch enrichment analysis
    """

    sb.set_theme(style="dark")
    plt.hist([actives_df[_AE.ACTIVES], inactives_df[_AE.INACTIVES]], color=["g", "b"],
             label=[_AE.ACTIVES, _AE.INACTIVES], bins=20, stacked=True, alpha=0.7)

    if os.path.isdir(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_ACTIVES]):
        file_name = csv_name.replace(".csv", "")
    else:
        file_name = parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_ACTIVES].rsplit("/", 1)[-1].replace(".csv", "")

    plt.title(f"{file_name} Enrichment Histogram")
    plt.xlabel("Docking Score")
    plt.ylabel("Absolute Counts")
    plt.legend()

    output_path = os.path.join(parameters[_AE.OUTPUT][_AE.OUTPUT_PATH], f"{file_name}_enrichment_histogram.png")
    plt.savefig(output_path, format="png", dpi=300)
    plt.figure()


def enrichment_boxplot(parameters: dict, actives_df: pd.DataFrame, inactives_df: pd.DataFrame, csv_name=None):
    """this function takes the parameters dictionary provided by the user containing all required parameters.
       The actives and inactives data is used to construct an enrichment box plot

    :param parameters: a dictionary containing all the required parameters provided by the user
    :param actives_df: a DataFrame containing the actives data
    :param inactives_df: a DataFrame containing the inactives data
    :param csv_name: the name of the current data file being analyzed. This parameter is only relevant during
                     batch enrichment analysis
    """

    # calculate statistic metrics
    actives_stats = actives_df.describe().drop(["std"]).round(decimals=2).to_string()
    inactives_stats = inactives_df.describe().drop(["std"]).round(decimals=2).to_string()

    sb.set_theme(style="dark")
    fig, ax = plt.subplots()
    ax.boxplot([actives_df[_AE.ACTIVES], inactives_df[_AE.INACTIVES]], widths=[0.5, 0.5],
               labels=[_AE.ACTIVES, _AE.INACTIVES])

    text_params = dict(boxstyle="round", facecolor="wheat", edgecolor="black", alpha=0.3)
    ax.text(0.01, 0.98, actives_stats, transform=ax.transAxes, fontsize=7, verticalalignment="top", bbox=text_params)
    ax.text(0.84, 0.98, inactives_stats, transform=ax.transAxes, fontsize=7, verticalalignment="top", bbox=text_params)

    if os.path.isdir(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_ACTIVES]):
        file_name = csv_name.replace(".csv", "")
    else:
        file_name = parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_ACTIVES].rsplit("/", 1)[-1].replace(".csv", "")

    plt.title(f"{file_name} Enrichment Boxplot")
    plt.xlabel("Ligand Activity")
    plt.ylabel("Docking Score")

    output_path = os.path.join(parameters[_AE.OUTPUT][_AE.OUTPUT_PATH], f"{file_name}_enrichment_boxplot.png")
    plt.savefig(output_path, format="png", dpi=300)
    plt.figure()


def pROC_curve(parameters:dict, TPR: list, FPR: list, rand_selection: list, FPR_AUC: list, pROC_AUC_values, csv_name=None) -> dict:
    """this function takes the parameters dictionary provided by the user containing all required parameters.
       The true positive rates (TPR), false positive rates (FPR), random selection rates (rand_selection)
       are used to construct a pROC curve. The false positive rates corresponding to when actives were
       recovered is given by FPR_AUC and is used to calculated the pROC AUC. A dictionary containing all
       the calculated pROC AUC values is returned


    :param parameters: a dictionary containing all the required parameters provided by the user
    :param TPR: true positive rates list
    :param FPR: false positive rates list
    :param rand_selection: random selection rates list
    :pROC_AUC_values: dictionary to store all calculated pROC AUC values
    :param csv_name: the name of the current data file being analyzed. This parameter is only relevant during
                     batch enrichment analysis
    """

    sb.set_theme(style="dark")
    fig, ax = plt.subplots()
    ax.plot(FPR, TPR, label="Docking Protocol")
    ax.plot(rand_selection, rand_selection, color="gray", label="Random")
    plt.legend()
    plt.ylim(0, 1.01)
    ax.margins(x=0.005)
    ax.set_xticks([0, 0.001, 0.010, 0.100, 1.000])
    ax.set_xscale("log")
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    if os.path.isdir(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_ACTIVES]):
        file_name = csv_name.replace(".csv", "")
    else:
        file_name = parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_ACTIVES].rsplit("/", 1)[-1].replace(".csv", "")

    plt.title(f"{file_name} pROC Curve")
    plt.xlabel("False Positive Rate (FPR)")
    plt.ylabel("True Positive Rate (TPR)")

    output_path = os.path.join(parameters[_AE.OUTPUT][_AE.OUTPUT_PATH], f"{file_name}_pROC_curve.png")
    plt.savefig(output_path, format="png", dpi=300)
    plt.figure()

    # calculate pROC AUC and store it in pROC_AUC_values
    pROC_AUC_values[file_name] = analysis.pROC_AUC(FPR_AUC)

    return pROC_AUC_values


def pROC_overlay_curve(parameters:dict, pROC_AUC_values, EF_values):
    """this function takes the parameters dictionary provided by the user containing all required parameters. An
       overlay pROC curve is generated and a dictionary containing all the calculated pROC AUC values are saved
       to the user specified output directory

    :param parameters: a dictionary containing all the required parameters provided by the user
    :param pROC_AUC_values: a dictionary to store all calculated pROC AUC values
    :param EF_values: a dictionary to store all calculated EF 5% values
    """

    sb.set_theme(style="dark")
    fig, ax = plt.subplots()

    csv_names = [csv for csv in os.listdir(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_ACTIVES]) if csv.endswith(".csv")]

    for idx, csv_name in enumerate(csv_names):
        # load the actives and inactives data
        actives_data = pd.read_csv(os.path.join(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_ACTIVES], csv_name))
        inactives_data = pd.read_csv(os.path.join(parameters[_AE.INPUT_ENRICHMENT_DATA][_AE.DATA_PATH_INACTIVES], csv_name))
        # generate pROC curve data points
        TPR, FPR, rand_selection, FPR_AUC = analysis.pROC_curve_datapoints(parameters, actives_data, inactives_data)

        file_name = csv_name.replace(".csv", "")

        if idx == 0:
            ax.plot(rand_selection, rand_selection, color="gray", label="Random")

        ax.plot(FPR, TPR, label=file_name)
        # calculate pROC AUC and store it in pROC_AUC_values
        pROC_AUC_values[file_name] = analysis.pROC_AUC(FPR_AUC)
        # calculate EF 5% and store it in EF_values
        EF_values = analysis.enrichment_factor(parameters, actives_data, inactives_data, EF_values, csv_name)

    plt.legend()
    plt.ylim(0, 1.01)
    ax.margins(x=0.005)
    ax.set_xticks([0, 0.001, 0.010, 0.100, 1.000])
    ax.set_xscale("log")
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    plt.title("Overlay pROC Curve")
    plt.xlabel("False Positive Rate (FPR)")
    plt.ylabel("True Positive Rate (TPR)")

    output_path = os.path.join(parameters[_AE.OUTPUT][_AE.OUTPUT_PATH], "overlay_pROC_curve.png")
    plt.savefig(output_path, format='png', dpi=300)
    plt.figure()

    # save the pROC AUC values in json format for readability
    json_output_path = os.path.join(parameters[_AE.OUTPUT][_AE.OUTPUT_PATH], "pROC_AUC_values.json")
    with open(json_output_path, "w+") as f:
        json.dump(pROC_AUC_values, f, indent=2)

    # save the EF 5% values in json format for readability
    json_output_path = os.path.join(parameters[_AE.OUTPUT][_AE.OUTPUT_PATH], "EF_5%_values.json")
    with open(json_output_path, "w+") as f:
        json.dump(EF_values, f, indent=2)

    sys.exit("Overlay pROC curves plot generated. All pROC AUC calculated. Exiting script now.")


def output_analysis(input_json: str):
    """this function takes a user provided json file containing all required analysis metrics and
       calls helper functions to perform statistical analyses and generate plots

    :param input_json: user provided json file containing all required metrics for analysis
    :raises FileNotFoundError: this error is raised if data_path or exp_data_path are neither a valid path
                               to a folder nor a file
    """

    # load the user input json into a dictionary
    with open(input_json, "r") as f:
        parameters = json.load(f)

    # call enrichment_analysis if user specifies for enrichment analysis
    enrichment_bool = parameters[_AE.PLOT_SETTINGS][_AE.ENRICHMENT_ANALYSIS]
    if enrichment_bool.lower() == "true":
        enrichment_bool = True
    elif enrichment_bool.lower() == "false":
        enrichment_bool = False
    else:
        sys.exit("'enrichment_analysis' value error: set to either 'True' or 'False'. For Enrichment Analysis, set to "
                 "'True'. For either Correlation or Thresholds Analysis, set to 'False'. Exiting script.")

    # check if output_path exists and if not, create it
    if not os.path.exists(parameters[_AE.OUTPUT][_AE.OUTPUT_PATH]):
        os.mkdir(parameters[_AE.OUTPUT][_AE.OUTPUT_PATH])

    if enrichment_bool:
        # initialize a dictionary to hold pROC AUC values
        pROC_AUC_values = {}
        # initialize a dictionary to hold EF 5% vaues
        EF_values = {}
        # add the random selection pROC AUC value
        pROC_AUC_values["Random"] = 0.434
        enrichment_analysis(parameters, pROC_AUC_values, EF_values)

    # at this point, the user wants either Correlation or Thresholds Analysis
    # first check if data_path is valid (either a folder or single csv file)
    if not os.path.isdir(parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_PATH]) and not os.path.isfile(parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_PATH]):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_PATH])
    # if data_path is a folder, ensure it is not empty and that it contains at least 1 csv file
    if os.path.isdir(parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_PATH]):
        if not os.listdir(parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_PATH]):
            sys.exit(parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_PATH] + " folder is empty. Please ensure your DockStream output csv files are added to the folder.")
        elif not any(file.endswith('.csv') for file in os.listdir(parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_PATH])):
            sys.exit(parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_PATH] + " contains no csv files. Please ensure your DockStream output csv files are added to the folder.")
    # at this point, data_path must be a file. Check that it is in csv format
    if os.path.isfile(parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_PATH]):
        if not parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_PATH].endswith(".csv"):
            sys.exit(parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_PATH] + " is not a csv file. Please ensure it is in csv format.")
    # check that exp_data_path is valid
    if not os.path.isfile(parameters[_AE.INPUT_EXP_DATA][_AE.EXP_DATA_PATH]):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), parameters[_AE.INPUT_EXP_DATA][_AE.EXP_DATA_PATH])
    # check that exp_data_path is a csv file
    if not parameters[_AE.INPUT_EXP_DATA][_AE.EXP_DATA_PATH].endswith(".csv"):
        sys.exit(parameters[_AE.INPUT_EXP_DATA][_AE.EXP_DATA_PATH] + " should be a csv file.")

    # go ahead with the rest of the script for Correlation/Thresholds Analysis
    # if the data path is a directory, extract all csv files
    if os.path.isdir(parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_PATH]):
        docking_results = [file for file in os.listdir(parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_PATH]) if file.endswith(".csv")]
    # the data path must be a single file then
    else:
        docking_results = [os.path.basename(os.path.normpath(parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_PATH]))]

    # initialize a dictionary to hold statistics metrics
    analysis_results = {}

    for csv_file in docking_results:

        # call correlation_analysis from analysis.py to determine the coefficient of determination
        # and Spearman and Kendall rank correlation coefficients
        parameters, analysis_results, comparison_df = analysis.correlation_analysis(parameters, csv_file, analysis_results)

        # check if the user specified docking and experimental thresholds are "N/A". In which case, perform
        # correlation analysis - generate scatter plots and compute correlation metrics only
        if parameters[_AE.INPUT_DOCKING_DATA][_AE.DATA_THRESHOLDS] == "N/A" and parameters[_AE.INPUT_EXP_DATA][_AE.EXP_THRESHOLDS] == "N/A":
            scat_plot(parameters, csv_file, comparison_df)
            continue

        # call scat_plot to create a scatter plot of the docking scores against the experimental data
        scat_plot(parameters, csv_file, comparison_df)

        # call binary_matrix to create a confusion matrix of the docking scores against the experimental data
        matt_coeffs = binary_matrix(parameters, csv_file, comparison_df)
        # loop through matt_coeffs and populate analysis_results with the MCC values
        for trial, mcc in matt_coeffs.items():
            analysis_results[trial] = mcc

        # call histogram to create histograms
        histogram(parameters, csv_file, comparison_df)

    # save the backend performance metrics in json format for readability
    json_output_path = os.path.join(parameters[_AE.OUTPUT][_AE.OUTPUT_PATH], "all_runs_stats_metrics.json")
    with open(json_output_path, 'w+') as f:
        json.dump(analysis_results, f, indent=2)


if __name__ == "__main__":
    # take user specified input parameters to run the analysis script
    parser = argparse.ArgumentParser(description="Implements entry point to DockStream output analysis.")
    parser.add_argument("-input_json", type=str, required=True, help="Path to user provided json file which contains all the paths/metrics for analysis.")
    args = parser.parse_args()

    output_analysis(args.input_json)
