from nolabs.application.jobs.use_cases import (
    DeleteJobFeature, GetJobMetadataFeature, GetJobsMetadataFeature,
    UpdateJobFeature)


class JobDependencies:
    @staticmethod
    def jobs_metadata() -> GetJobsMetadataFeature:
        return GetJobsMetadataFeature()

    @staticmethod
    def job_metadata() -> GetJobMetadataFeature:
        return GetJobMetadataFeature()

    @staticmethod
    def delete_job():
        return DeleteJobFeature()

    @staticmethod
    def update_job():
        return UpdateJobFeature()
