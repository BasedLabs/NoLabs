from reinvent_microservice import Configuration, ApiClient, DefaultApi

from nolabs.api_models.small_molecules_design import GetExperimentResponse, ExperimentPropertiesResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.infrastructure.settings import Settings
from nolabs.modules.small_molecules_design.services.file_management import FileManagement


class GetExperimentFeature:
    def __init__(self, settings: Settings, file_management: FileManagement):
        self._settings = settings
        self._fm = file_management

    async def handle(self, experiment_id: str) -> GetExperimentResponse | None:
        if not self._fm.experiment_exists(experiment_id):
            return None

        configuration = Configuration(
            host=self._settings.reinvent_host,
        )

        experiment = self._fm.get_metadata(ExperimentId(experiment_id))

        params = self._fm.get_params(ExperimentId(experiment_id))
        receptor = self._fm.get_pdb(ExperimentId(experiment_id))

        if not receptor:
            properties = ExperimentPropertiesResponse(
                pdb_file=None,
                pdb_file_name=None,
                center_x=0.0,
                center_y=0.0,
                center_z=0.0,
                size_x=0.0,
                size_y=0.0,
                size_z=0.0,
                batch_size=128,
                minscore=0.4,
                epochs=50
            )
        else:
            properties = ExperimentPropertiesResponse(
                pdb_file=receptor.read_string(),
                pdb_file_name=receptor.name,
                center_x=params.center_x,
                center_y=params.center_y,
                center_z=params.center_z,
                size_x=params.size_x,
                size_y=params.size_y,
                size_z=params.size_z,
                batch_size=params.batch_size,
                minscore=params.minscore,
                epochs=params.epochs,
            )

        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)

            job_result = api_instance.get_job_jobs_job_id_get(experiment_id)

            running = job_result.running if job_result else False
            learning_completed = job_result.learning_completed if job_result else False

            return GetExperimentResponse(
                experiment_id=experiment_id,
                experiment_name=experiment.name.value,
                created_at=experiment.date,
                running=running,
                learning_completed=learning_completed,
                properties=properties
            )
