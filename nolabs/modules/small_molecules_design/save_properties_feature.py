from reinvent_microservice import Configuration, ApiClient, DefaultApi

from nolabs.api_models.small_molecules_design import ExperimentPropertiesRequest
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.small_molecules_design.services.file_management import FileManagement
from nolabs.infrastructure.settings import Settings


class SavePropertiesFeature:
    def __init__(self, file_management: FileManagement, settings: Settings):
        self._settings = settings
        self._fm = file_management

    async def handle(self, experiment_id: str, request: ExperimentPropertiesRequest):
        configuration = Configuration(
            host=self._settings.reinvent_host,
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)

            api_instance.save_params_jobs_job_id_params_post(
                job_id=request.job_id,
                pdb_file=request.pdb_file,
                center_x=request.center_x,
                center_y=request.center_y,
                center_z=request.center_z,
                size_x=request.size_x,
                size_y=request.size_y,
                size_z=request.size_z,
                batch_size=request.batch_size,
                minscore=request.minscore,
                epochs=request.epochs
            )

        self._fm.set_params(ExperimentId(experiment_id), params=request)
