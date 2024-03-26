import datetime
import os
import subprocess
import sys
import tomllib

import pydantic
import toml

import uuid
from typing import List, Any

from fastapi import UploadFile

from leaf import FileObject, DirectoryObject
from starlette.responses import FileResponse

from microservice.exceptions import ReinventException, ErrorCode
from microservice.api_models import RunFineTuningJobRequest, FineTuningJobResponse


@pydantic.dataclasses.dataclass
class Process:
    id: str
    running: bool
    started_at: datetime.datetime
    progress: float
    pdbqt_filename: str
    epochs: int
    directory: DirectoryObject
    errors_file: FileObject
    output_file: FileObject
    process: Any


processes: List[Process] = []


def get_process(id) -> Process | None:
    p = [p for p in processes if p if p.id == id]
    if p:
        return p[0]
    return None


class FineTuning:
    def __init__(self):
        self.fs = DirectoryObject(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'finetuning'))

    def get_job(self, id: str) -> FineTuningJobResponse | None:
        directory = self.fs.directories.first_or_default(lambda o: o.name == id)
        process = get_process(id)
        if not process:
            return None
        pdbqt_file = directory.files.first_or_default(lambda o: o.name == process.pdbqt_filename,
                                                      recursive=True)
        pdbqt_content = pdbqt_file.read_string()

        return FineTuningJobResponse(
            id=id,
            running=process.running,
            started_at=process.started_at,
            progress=process.progress,
            pdbqt_content=pdbqt_content,
            pdbqt_filename=process.pdbqt_filename,
            epochs=process.epochs,
            errors=process.errors_file.read_string()
        )

    def calculate_progress(self, metadata, current_stage) -> float:
        return float(current_stage) / float(metadata['epochs'])

    def all_jobs(self) -> List[FineTuningJobResponse]:
        result = []
        for run_dir in self.fs.children:
            job = self.get_job(run_dir.name)
            if job:
                result.append(job)
        return result

    async def prepare_pdbqt(self, pdb: UploadFile) -> FileResponse:
        tmp = self.fs.add_directory('tmp')
        pdb_file = tmp.add_file(str(uuid.uuid4()) + '.pdb')
        pdb_file.write_bytes(await pdb.read())
        pdbqt = self.fs.add_directory('tmp').add_file(str(uuid.uuid4()) + '.pdbqt')
        prepare_pdbqt = self.fs.files.first_or_default(lambda o: o.name == 'prepare_pdbqt.sh')
        subprocess.run([prepare_pdbqt.full_path, pdb_file.full_path, pdbqt.full_path])
        return FileResponse(pdbqt.full_path)

    def run_reinforcement_learning(self, RL_toml_file: FileObject, log_file: FileObject, error_file: FileObject) -> \
            subprocess.Popen[bytes]:
        """Returns PID of the process"""
        run_reinforcement_learning_shell = self.fs.files.first_or_default(
            lambda o: o.name == 'start_reinforcement_learning.sh')
        print('FULL PATH', [run_reinforcement_learning_shell.full_path,
                            RL_toml_file.full_path,
                            log_file.full_path,
                            error_file.full_path,
                            ])
        return subprocess.Popen([run_reinforcement_learning_shell.full_path,
                                 RL_toml_file.full_path,
                                 # log_file.full_path,
                                 # error_file.full_path,
                                 ], shell=False, stdout=sys.stdout, stderr=sys.stderr)

    def prepare_dockstream_config(self,
                                  request: RunFineTuningJobRequest,
                                  pdbqt: FileObject,
                                  directory: DirectoryObject) -> FileObject:
        # Configure dockstream
        dockstream_config = directory.files.first_or_default(lambda o: o.name == 'dockstream.json')

        dockstream_config_json = dockstream_config.read_json()
        dockstream_config_json['docking']['docking_runs'][0]['parameters']['parallelization'][
            'number_cores'] = os.cpu_count() * 2  # type: ignore
        dockstream_config_json['docking']['docking_runs'][0]['parameters']['receptor_pdbqt_path'].append(
            pdbqt.full_path)
        dockstream_config_json['docking']['docking_runs'][0]['parameters']['number_poses'] = 2

        dockstream_config_json['docking']['docking_runs'][0]['parameters']['search_space'][
            '--center_x'] = request.center_x
        dockstream_config_json['docking']['docking_runs'][0]['parameters']['search_space'][
            '--center_y'] = request.center_y
        dockstream_config_json['docking']['docking_runs'][0]['parameters']['search_space'][
            '--center_z'] = request.center_z
        dockstream_config_json['docking']['docking_runs'][0]['parameters']['search_space']['--size_x'] = request.size_x
        dockstream_config_json['docking']['docking_runs'][0]['parameters']['search_space']['--size_y'] = request.size_y
        dockstream_config_json['docking']['docking_runs'][0]['parameters']['search_space']['--size_z'] = request.size_z

        poses = directory.add_file('poses.sdf')
        scores = directory.add_file('scores.sdf')
        dockstream_config_json['docking']['docking_runs'][0]['output']['poses']['poses_path'] = poses.full_path
        dockstream_config_json['docking']['docking_runs'][0]['output']['scores']['scores_path'] = scores.full_path

        print('CONFIG', dockstream_config_json)

        dockstream_config.write_json(dockstream_config_json)

        return dockstream_config

    def prepare_toml(self, epochs: int, dockstream_config: FileObject, directory: DirectoryObject) -> FileObject:
        # Configure learning
        reinforcement_learning_config = directory.files.first_or_default(lambda o: o.name == 'RL.toml')
        t = reinforcement_learning_config.read(tomllib.load, mode='rb')
        t['output_csv'] = directory.add_file('rl_direct.csv').full_path
        t['stage'][0]['chkpt_file'] = directory.add_file('rl_direct.chkpt').full_path
        t['stage'][0]['max_steps'] = epochs
        t['stage'][0]['scoring']['component'][0]['DockStream']['endpoint'][0] \
            ['params']['configuration_path'] = dockstream_config.full_path

        reinforcement_learning_config.write_string(toml.dumps(t))

        return reinforcement_learning_config

    async def run(self, pdbqt_file: UploadFile, request: RunFineTuningJobRequest):
        if not pdbqt_file:
            raise ReinventException(ErrorCode.PDBQT_NOT_PROVIDED)

        id = str(uuid.uuid4())
        directory = self.fs.copy(self.fs.add_directory(id))
        directory.directories.first_or_default(lambda o: o.name == id).delete()
        pdbqt = directory.add_file(pdbqt_file.filename)
        pdbqt.write_bytes(await pdbqt_file.read())

        dockstream_config_json = self.prepare_dockstream_config(request, pdbqt, directory)
        reinforcement_learning_config = self.prepare_toml(request.epochs, dockstream_config_json, directory)

        print('TOML', reinforcement_learning_config.read_string())

        log_file = directory.add_file('output.log')
        error_file = directory.add_file('error.log')
        popen = self.run_reinforcement_learning(reinforcement_learning_config, log_file, error_file)

        processes.append(Process(
            id=id,
            running=True,
            started_at=datetime.datetime.now(datetime.UTC),
            progress=0.0,
            pdbqt_filename=pdbqt_file.filename,  # type: ignore
            epochs=request.epochs,
            directory=directory,
            errors_file=error_file,
            output_file=log_file,
            process=popen))

        return self.get_job(id)
