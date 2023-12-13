import os
import json

class ProgressTracker:
    def __init__(self, experiment_dir, protein_ids, model_tasks):
        self.experiment_dir = experiment_dir
        self.protein_ids = protein_ids
        self.model_tasks = model_tasks
        self.progress_file = os.path.join(experiment_dir, "progress.json")
        self._initialize_progress_file()

    def _initialize_progress_file(self):
        if not os.path.exists(self.progress_file):
            initial_progress = {
                "experiment_progress": 0,
                "total_sequences": len(self.protein_ids),
                "processed_sequences": 0
            }
            with open(self.progress_file, "w") as f:
                json.dump(initial_progress, f)

        # Initialize or load individual protein progress files
        for protein_id in self.protein_ids:
            protein_dir = os.path.join(self.experiment_dir, protein_id)
            os.makedirs(protein_dir, exist_ok=True)
            protein_progress_file = os.path.join(protein_dir, "progress.json")
            if not os.path.exists(protein_progress_file):
                protein_progress = {
                    "completed_tasks": [],
                    "all_tasks": self.model_tasks,
                    "progress": 0
                }
                with open(protein_progress_file, "w") as f:
                    json.dump(protein_progress, f)

    def update_protein_progress(self, protein_id, model_task):
        protein_progress_file = os.path.join(self.experiment_dir, protein_id, "progress.json")
        with open(protein_progress_file, "r+") as f:
            protein_progress = json.load(f)
            if model_task not in protein_progress["completed_tasks"]:
                protein_progress["completed_tasks"].append(model_task)
                protein_progress["progress"] = 100*len(protein_progress['completed_tasks'])/len(protein_progress['all_tasks'])
                f.seek(0)
                f.truncate()
                json.dump(protein_progress, f)

    def update_experiment_progress(self):
        with open(self.progress_file, "r+") as f:
            progress_data = json.load(f)
            progress_data["processed_sequences"] += 1
            progress_data["experiment_progress"] = (progress_data["processed_sequences"] / len(self.protein_ids)) * 100
            f.seek(0)
            f.truncate()
            json.dump(progress_data, f)

def get_progress(experiment_dir):
    progress_file = os.path.join(experiment_dir, "progress.json")
    with open(progress_file, "r") as f:
        progress_data = json.load(f)
    return {"progress": progress_data["experiment_progress"]}

def get_protein_progress(experiment_dir, protein_id):
    protein_progress_file = os.path.join(experiment_dir, protein_id, "progress.json")
    with open(protein_progress_file, "r") as f:
        protein_progress = json.load(f)
    return protein_progress
