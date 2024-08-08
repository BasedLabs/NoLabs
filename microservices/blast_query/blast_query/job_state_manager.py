class JobState:
    def __init__(self):
        self.is_running = False


class JobStateManager:
    def __init__(self):
        self.jobs = {}  # Maps job_id to JobState

    def start_job(self, job_id: str):
        if job_id not in self.jobs:
            self.jobs[job_id] = JobState()
        self.jobs[job_id].is_running = True

    def finish_job(self, job_id: str):
        if job_id in self.jobs:
            del self.jobs[job_id]

    def is_job_running(self, job_id: str) -> bool:
        return self.jobs.get(job_id, JobState()).is_running

    def get_running_jobs(self):
        return [job_id for job_id, state in self.jobs.items() if state.is_running]


# Global instance
job_state_manager = JobStateManager()
