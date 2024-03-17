import datetime
import os
import subprocess

import leaf as lf
import uuid
from typing import List, Any, Dict

from fastapi import UploadFile

from leaf.exception import ErrorCodes
from leaf import FileObject, DirectoryObject
from leaf.filesystem import ObjectType
from microservice.exceptions import ReinventException, ErrorCode
from microservice.api_models import RunFineTuningJobRequest, FineTuningJobResponse


class FineTuning:
    def __init__(self):
        self.fs = lf.DirectoryObject(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'finetuning'))

    def get(self, id: str) -> FineTuningJobResponse | None:
        metadata_obj = (self.fs
                        .files
                        .first_or_default(lambda o: o.parent.name == id and o.name == 'metadata.json',
                                             recursive=True))

        if not metadata_obj:
            return

        metadata: Dict[str, Any] = metadata_obj.read_json()
        pdb_filename = metadata['pdb_filename']
        pdb_file: lf.FileObject | None = self.fs.files.first_or_default(lambda o: o.name == pdb_filename, recursive=True)
        pdb_content = pdb_file.read_string()

        return FineTuningJobResponse(
            id=id,
            name=metadata['name'],
            running=metadata['running'],
            started_at=metadata['started_at'],
            finished_at=metadata['finished_at'],
            progress=metadata['progress'],
            pdb_content=pdb_content,
            pdb_filename=pdb_filename,
            errors=metadata['errors'],
            epochs=metadata['epochs']
        )

    def refresh_progress(self):
        job_directory: DirectoryObject
        for job_directory in self.fs.where(lambda o: o.type == ObjectType.DIRECTORY):
            metadata_file: FileObject = job_directory.files.first_or_default(lambda o: o.name == 'metadata.json')
            metadata = metadata_file.read_json()
            pid = metadata['pid']
            log_file: FileObject = job_directory.files.first_or_default(lambda o: o.name == f'{pid}.log')
            errors_file: FileObject = job_directory.files.first_or_default(lambda o: o.name == f'{pid}.error')
            metadata['errors'] = errors_file.read_string()


    def all(self) -> List[FineTuningJobResponse]:
        result = []
        for run_dir in self.fs.children:
            job = self.get(run_dir.name)
            if job:
                result.append(job)
        return result

    def prepare_pdbqt(self, pdb_file: FileObject, pdbqt: FileObject):
        prepare_pdbqt = self.fs.files.first_or_default(lambda o: o.name == 'prepare_pdbqt.sh')
        subprocess.run([prepare_pdbqt.full_path, pdb_file.full_path, pdbqt.full_path])

    async def run(self, pdb_file: UploadFile, request: RunFineTuningJobRequest):
        if not pdb_file:
            raise ReinventException(ErrorCode.PDB_NOT_PROVIDED)

        id = str(uuid.uuid4())
        directory = self.fs.add_directory(id)
        pdb = directory.add_file(pdb_file.filename)
        pdb.write_bytes(await pdb_file.read())

        # Configure dockstream
        dockstream_config = self.fs.files.first_or_default(lambda o: o.name == 'dockstream.json')
        pdbqt = directory.add_file(pdb.name + 'qt')
        self.prepare_pdbqt(pdb, pdbqt)

        dockstream_config_json = dockstream_config.read_json()
        dockstream_config_json['docking_runs'][0]['parameters']['parallelization']['number_cores'] = os.cpu_count()
        dockstream_config_json['docking_runs'][0]['parameters']['receptor_pdbqt_path'].append(pdbqt.full_path)
        dockstream_config_json['docking_runs'][0]['parameters']['number_poses'] = 2

        dockstream_config_json['docking_runs'][0]['parameters']['search_space']['--center_x'] = request.center_x
        dockstream_config_json['docking_runs'][0]['parameters']['search_space']['--center_y'] = request.center_y
        dockstream_config_json['docking_runs'][0]['parameters']['search_space']['--center_z'] = request.center_z
        dockstream_config_json['docking_runs'][0]['parameters']['search_space']['--size_x'] = request.size_x
        dockstream_config_json['docking_runs'][0]['parameters']['search_space']['--size_y'] = request.size_y
        dockstream_config_json['docking_runs'][0]['parameters']['search_space']['--size_z'] = request.size_z

        poses = directory.add_file('poses.sdf')
        scores = directory.add_file('scores.sdf')
        dockstream_config_json['docking_runs'][1]['poses']['poses_path'] = poses.full_path
        dockstream_config_json['docking_runs'][1]['scores']['scores_path'] = scores.full_path

        # Configure learning
        reinforcement_learning_config = self.fs.files.first_or_default(lambda o: o.name == 'RL.toml')

        directory.add_file('metadata.json').write_json(
            {
                'name': request.name,
                'running': False,
                'started_at': datetime.datetime.now(datetime.UTC),
                'finished_at': None,
                'progress': 0.0,
                'pdb_filename': pdb_file.filename,
                'errors': '',
                'epochs': request.epochs
            }
        )

