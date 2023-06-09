import os
import csv
from datetime import date
from typing import List

from .model import BaseModel


class Pipeline:
    def __init__(self, models: List[BaseModel] = []):
        """
        models: pass the list of loaded models which inherit BaseModel class
        """
        self.models = models

    def get_model_names(self):
        return [model.model_name for model in self.models]

    def get_model_by_task(self, model_task: str):
        """
        Gets the first model of the required type from the pipeline (TODO: get all the models of correct type)
        """
        for model in self.models:
            if model.model_task == model_task:
                return model

    def add_model(self, model: BaseModel):
        """
        model: loaded model which inherits the BaseModel class
        """
        self.models.append(model)

    def predict(self, sequence: str):
        """
        Make predictions for one sequence
        """
        output = {}
        for model in self.models:
            output[model.model_name] = model.predict(sequence)
        return output

    def predict_multiple(self, sequences: List[str], save_to_csv=False, custom_filename=""):
        """
        Make predictions for multiple sequences
        """
        results = []
        for sequence in sequences:
            output = self.predict(sequence)
            results.append((sequence, output))

        if save_to_csv:
            # Currently will save into results/ directory
            save_directory = os.path.dirname(os.path.abspath(__file__)) + "/results/"
            check_folder_exists(save_directory)
            if custom_filename:
                final_path = save_directory + custom_filename + ".csv"
                self.save_results_to_csv(results, final_path)
            else:
                final_path = save_directory + str(date.today()) + ".csv"
                self.save_results_to_csv(results, final_path)
        return results

    def save_results_to_csv(self, results, filename):
        if not results:
            print(f"Your predictions are empty")
            return

        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['Sequence'] + self.get_model_names())

            for result in results:
                row = [result[0]]
                outputs = result[1]

                for key in outputs.keys():
                    row.append(outputs[key])

                writer.writerow(row)


def check_folder_exists(filename):
    directory = os.path.dirname(filename)

    if os.path.exists(directory):
        print(f"The folder for '{filename}' does not exist. Creating one")
        os.makedirs(directory)
