from typing import TypeVar, Generic

from nolabs.api_models.amino_acid.common_models import RunAminoAcidRequest
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.modules.amino_acid.file_management_base import AminoAcidFileManagementBase


TFileManagement = TypeVar('TFileManagement', bound=AminoAcidFileManagementBase, contravariant=True)

class RunAminoAcidInferenceFeature(Generic[TFileManagement]):
    def __init__(self, file_management: TFileManagement):
        self._file_management = file_management

    async def _setup_experiment(self,
                                experiment_id: ExperimentId,
                                request: RunAminoAcidRequest):
        self._file_management.ensure_experiment_folder_exists(experiment_id=experiment_id)
        self._file_management.cleanup_experiment(experiment_id=experiment_id)
        await self._file_management.set_metadata(experiment_id=experiment_id,
                                                 experiment_name=ExperimentName(request.experiment_name))
        await self._file_management.set_properties(experiment_id=experiment_id, request=request)
