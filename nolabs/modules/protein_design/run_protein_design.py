import protein_design_microservice as microservice

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.api_models.protein_design import RunProteinDesignRequest, RunProteinDesignResponse
from nolabs.infrastructure.settings import Settings
from nolabs.modules.protein_design.services.file_management import FileManagement


class RunProteinDesignFeature:
    def __init__(self,
                 file_management: FileManagement,
                 settings: Settings):
        self._file_management = file_management
        self._settings = settings

    async def handle(self,
                     request: RunProteinDesignRequest) -> RunProteinDesignResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        experiment_name = ExperimentName(request.experiment_name)
        pdb_content = await request.pdb_file.read()

        self._file_management.ensure_experiment_folder_exists(experiment_id=experiment_id)
        self._file_management.cleanup_experiment(experiment_id=experiment_id)
        await self._file_management.set_metadata(experiment_id=experiment_id, experiment_name=experiment_name)
        await self._file_management.set_properties(experiment_id=experiment_id, request=request)

        configuration = microservice.Configuration(
            host=self._settings.protein_design_host
        )
        with microservice.ApiClient(configuration=configuration) as client:
            api_instance = microservice.DefaultApi(client)
            response = api_instance.run_rfdiffusion_endpoint_run_rfdiffusion_post(
                run_rfdiffusion_request=microservice.RunRfdiffusionRequest(
                    pdb_content=pdb_content.decode('utf-8'),
                    hotspots=request.hotspots,
                    contig=request.contig,
                    timesteps=request.timesteps,
                    number_of_designs=request.number_of_designs
                )
            )

            if response.errors and not response.pdbs_content:
                raise NoLabsException(response.errors, ErrorCodes.protein_design_run_error)

            await self._file_management.set_result(experiment_id=experiment_id, pdbs_content=response.pdbs_content,
                                                   request=request)
            return RunProteinDesignResponse(
                experiment_id=experiment_id.value,
                experiment_name=experiment_name.value,
                pdb_files=response.pdbs_content
            )