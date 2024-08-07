from nolabs.application.use_cases.jobs.use_cases import GetJobsMetadataFeature, DeleteJobFeature, \
    UpdateJobFeature, GetJobMetadataFeature


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