import os
import argparse
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy.stats import spearmanr
from scipy.stats import kendalltau
from matplotlib import pyplot as plt
import seaborn as sns


def to_scatter_plot(file_name: str, docking_scores: np.ndarray, experimental_scores: np.ndarray):

    sns.set()
    plt.scatter(docking_scores, experimental_scores, color='g')
    plt.title(file_name + ' Scatter Plot')
    plt.xlabel('Docking Score (-kcal/mol)')
    plt.ylabel('Experimental Binding (-kcal/mol)')
    plt.savefig(file_name.replace('csv', ''), dpi=300)
    plt.figure()


def output_analysis(data_dir: str, exp_data_path: str, output_dir: str) -> dict:

    docking_results = [file for file in os.listdir(data_dir) if file.endswith('.csv')]
    analysis_results = {}

    for file_name in docking_results:
        data_path = os.path.join(data_dir, file_name)
        # load the docking results
        docking_data = pd.read_csv(data_path)
        # load the experimental comparison data
        experimental_data = pd.read_csv(exp_data_path)
        # merge the docking results and experimental binding data based on their 'smiles' identifiers.
        # this is to ensure that data pertaining to the correct 'smiles' ligands are compared irrespective
        # of the order in which the user provides the experimental data
        comparison_df = docking_data.merge(experimental_data, on='smiles')

        # extract the DockStream docking scores and experimental potency parameter columns
        docking_scores = comparison_df['score'].to_numpy().reshape(-1,1)
        experimental_data = comparison_df['exp_binding'].to_numpy()
        # fit a linear model using least squares
        model = LinearRegression().fit(docking_scores, experimental_data)
        coeff_determination = model.score(docking_scores, experimental_data)
        # perform Spearman correlation analysis
        spearman_correlation = spearmanr(docking_scores, experimental_data)
        # perform Kendall correlation analysis
        kendall_correlation = kendalltau(docking_scores, experimental_data)
        # store the linear model and spearman quantities in a dictionary
        analysis_results[file_name] = {'coefficient of determination': coeff_determination, 'Spearman correlation':
                                        spearman_correlation, 'Kendall correlation': kendall_correlation}

        # call to_scatter_plot to create a scatter plot of docking_scores against experimental_data
        output_path = os.path.join(output_dir, file_name)
        to_scatter_plot(output_path, docking_scores, experimental_data)

    return analysis_results


if __name__ == '__main__':
    # take user specified input parameters to run the analysis script
    parser = argparse.ArgumentParser(description='Implements entry point to results analysis.')
    parser.add_argument('-data_dir', type=str, required=True, help='The path to the output csv files to be analyzed.')
    parser.add_argument('-exp_data_path', type=str, required=True, help='The path to the experimental binding data of your ligand library. This will be used for regression analysis.')
    parser.add_argument('-output_dir', type=str, required=True, help='The desired output folder path to store the analyzed results and plots.')
    #parser.add_argument('-regression_threshold', type=int, default=0.75, help='The desired regression threshold to classify good/poor backend performances based on your receptor-ligand system. The default value if not specified is 0.75.')
    args = parser.parse_args()

    analysis_results = output_analysis(args.data_dir, args.exp_data_path, args.output_dir)
