import os
import json
import errno
import sys
import argparse

from dockstream.utils.execute_external.execute import Executor
from dockstream.utils import files_paths
from dockstream.utils.enums.docking_enum import DockingConfigurationEnum

_DC = DockingConfigurationEnum()


def run_script(input_path: str) -> dict:
    """this method takes an input path to either a folder containing DockStream json files or a single json file and
       returns a dictionary whose keys are the json names and the corresponding values are the paths to the json
       file. The dictionary will be looped later to run DockStream

    :param input_path: path to either a folder of json files or a single json file
    :raises FileNotFoundError: this error is raised if input_path is neither a folder nor a file
    :return: dictionary, keys are the DockStream json names and values are the paths to them
    """

    # first check if input_path is valid (either a folder containing json files or a single json file)
    if not os.path.isdir(input_path) and not os.path.isfile(input_path):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), input_path)
    # if input_path is a folder, ensure it is not empty and that it contains at least 1 json file
    if os.path.isdir(input_path):
        if not os.listdir(input_path):
            sys.exit(input_path + ' folder is empty. Please ensure your DockStream json files are added to the folder.')
        elif not any(file.endswith('.json') for file in os.listdir(input_path)):
            sys.exit(input_path + ' contains no json files. Please ensure your DockStream json files are added to the folder.')
    # at this point, the path must be a file. Check that it is in json format
    if os.path.isfile(input_path):
        if not input_path.endswith('.json'):
            sys.exit(input_path + ' is not a json file. Please ensure it is in json format.')

    # initialize a dictionary to hold all DockStream runs
    batch_runs = {}
    # loop through all json files and update the paths if input_path if a directory
    if os.path.isdir(input_path):
        all_runs = [file for file in os.listdir(input_path) if file.endswith('.json')]
        for json in all_runs:
            batch_runs[json.replace('.json', '')] = os.path.join(input_path, json)

    # at this point, input path must be a single json file
    else:
        json_name = os.path.basename(os.path.normpath(input_path)).replace('.json', '')
        batch_runs[json_name] = input_path

    return batch_runs


if __name__ == '__main__':
    # take user specified input parameters to run the benchmarking script
    parser = argparse.ArgumentParser(description='Facilitates batch DockStream execution.')
    parser.add_argument('-input_path', type=str, required=True, help='The path to either a folder of DockStream json files or a single json file.')
    args = parser.parse_args()

    batch_runs = run_script(args.input_path)

    executor = Executor()
    # initialize a dictionary to store the names of all runs that did not enforce "best_per_ligand"
    non_bpl_runs = {}
    # loop through all user json files and run DockStream
    for trial_name, json_path in batch_runs.items():
        # check if the current DockStream run has "best_per_ligand" enforced
        with open(json_path, "r") as f:
            parameters = json.load(f)
            # in case output mode was not specified in the configuration json
            try:
                for docking_run in parameters[_DC.DOCKING][_DC.DOCKING_RUNS]:
                    output_mode = docking_run[_DC.OUTPUT][_DC.OUTPUT_SCORES][_DC.OUTPUT_MODE]
                    if output_mode != _DC.OUTPUT_MODE_BESTPERLIGAND:
                        non_bpl_runs[trial_name] = output_mode
                        break
            except:
                pass

        print(f'Running {trial_name}')

        result = executor.execute(command=sys.executable, arguments=[files_paths.attach_root_path('docker.py'),
                                  '-conf', json_path, '-debug'], check=False)

        print(result)
        # print out error messages (if applicable) for the current DockStream run
        if result.returncode != 0:
            print(f'There was an error with {trial_name} DockStream run.')
            print(result.stdout)
            print(result.stderr)

    if bool(non_bpl_runs):
        # print the names of the runs which did not enforce "best_per_ligand"
        print(f"List of runs which did not have 'best_per_ligand' specified. These runs cannot be "
              f"passed into the analysis script. {non_bpl_runs}")
