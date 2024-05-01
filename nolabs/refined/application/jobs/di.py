from nolabs.refined.application.jobs.use_cases import GetJobsMetadataFeature, DeleteJobFeature, \
    UpdateJobFeature


class JobDependencies:
    @staticmethod
    def job_metadata() -> GetJobsMetadataFeature:
        return GetJobsMetadataFeature()

    @staticmethod
    def delete_job():
        return DeleteJobFeature()

    @staticmethod
    def update_job():
        return UpdateJobFeature()