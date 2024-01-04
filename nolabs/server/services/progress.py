import os
import json

class ProgressTracker:
    def __init__(self, target_dir, tasks):
        self.target_dir = target_dir
        self.tasks = tasks
        self.progress_file = os.path.join(target_dir, "progress.json")
        self._initialize_progress_file()

    def _initialize_progress_file(self):
        if not os.path.exists(self.progress_file):
            initial_progress = {
                "progress": 0,
                "tasks": self.tasks,
                "completed_tasks": []
            }
            directory = os.path.dirname(self.progress_file)
            # Create the directory (and any intermediate directories) if it doesn't exist
            if not os.path.exists(directory):
                os.makedirs(directory)
            with open(self.progress_file, "w") as f:
                json.dump(initial_progress, f)

    def update_progress(self, completed_task):
        with open(self.progress_file, "r+") as f:
            progress = json.load(f)
            if completed_task not in progress["completed_tasks"]:
                progress["completed_tasks"].append(completed_task)
                progress["progress"] = 100*len(progress['completed_tasks'])/len(progress['tasks'])
                f.seek(0)
                f.truncate()
                json.dump(progress, f)

def get_progress(target_dir):
    progress_file = os.path.join(target_dir, "progress.json")
    try: 
        with open(progress_file, "r") as f:
            progress_data = json.load(f)
        return {"progress": progress_data["progress"]}
    except:
        return {"progress": 0}