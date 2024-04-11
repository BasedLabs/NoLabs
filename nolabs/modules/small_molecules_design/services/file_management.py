from typing import Optional, List

import pydantic
from leaf import DirectoryObject, FileObject

from nolabs.api_models.small_molecules_design import ExperimentPropertiesRequest
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.file_management_base import ExperimentsFileManagementBase
from nolabs.infrastructure.settings import Settings


@pydantic.dataclasses.dataclass
class ReinventParams:
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float
    batch_size: float
    minscore: float
    epochs: float


class FileManagement(ExperimentsFileManagementBase):
    def __init__(self, settings: Settings):
        super().__init__(settings.small_molecules_experiments_folder, settings.small_molecules_metadata_file_name)

    def get_jobs_ids(self, experiment_id: str) -> List[str]:
        ef_path = self.experiment_folder(ExperimentId(experiment_id))
        experiment_folder = DirectoryObject(ef_path)
        return [x.name for x in experiment_folder.directories]

    def get_pdb(self, experiment_id: ExperimentId) -> FileObject | None:
        ef_path = self.experiment_folder(experiment_id)
        experiment_folder = DirectoryObject(ef_path)

        return experiment_folder.files.first_or_default(lambda o: o.name == 'target.pdb')

    def save_pdb(self, experiment_id: ExperimentId, pdb: bytes):
        ef_path = self.experiment_folder(experiment_id)
        experiment_folder = DirectoryObject(ef_path)
        experiment_folder.add_file('target.pdb').write_bytes(pdb)

    def set_params(self, experiment_id: ExperimentId, params: ExperimentPropertiesRequest):
        ef_path = self.experiment_folder(experiment_id)
        experiment_folder = DirectoryObject(ef_path)
        experiment_folder.add_file('params.json').write_json({
            'center_x': params.center_x,
            'center_y': params.center_y,
            'center_z': params.center_z,
            'size_x': params.size_x,
            'size_y': params.size_y,
            'size_z': params.size_z,
            'batch_size': params.batch_size,
            'minscore': params.minscore,
            'epochs': params.epochs
        })

    def get_params(self, experiment_id: ExperimentId) -> Optional[ReinventParams]:
        ef_path = self.experiment_folder(experiment_id)
        experiment_folder = DirectoryObject(ef_path)
        params = experiment_folder.files.first_or_default(lambda o: o.name == 'params.json')
        if not params:
            return None
        j = params.read_json()

        return ReinventParams(
            center_x=j['center_x'],
            center_y=j['center_y'],
            center_z=j['center_z'],
            size_x=j['size_x'],
            size_y=j['size_y'],
            size_z=j['size_z'],
            batch_size=j['batch_size'],
            minscore=j['minscore'],
            epochs=j['epochs']
        )