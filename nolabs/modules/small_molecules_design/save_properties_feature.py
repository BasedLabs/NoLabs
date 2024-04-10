from reinvent_microservice import Configuration, ApiClient, ReinventApi

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

        pdb_file = await request.pdb_file.read()
        self._fm.save_pdb(ExperimentId(experiment_id), pdb_file)
        metadata = self._fm.get_metadata(ExperimentId(experiment_id))
        tmp_pdb = self._fm.get_pdb(ExperimentId(experiment_id))

        with ApiClient(configuration=configuration) as client:
            api_instance = ReinventApi(client)

            api_instance.save_params_api_reinvent_config_id_params_post(
                config_id=experiment_id,
                name=metadata.name.value,
                pdb_file=tmp_pdb.full_path,
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
