import dataclasses

import datetime
import os
import subprocess
import sys
import csv

import toml

import uuid
from typing import List, Any, Optional

from fastapi import UploadFile

from leaf import FileObject, DirectoryObject
from starlette.responses import FileResponse

from microservice.api_models import JobResponse, ParamsRequest, LogsResponse, ParamsResponse, Smiles, SmilesResponse
from microservice.exceptions import ReinventException, ErrorCode


@dataclasses.dataclass
class Process:
    id: str
    name: str
    running: bool
    learning_completed: bool
    created_at: datetime.datetime
    pdbqt_content: str
    pdbqt_filename: str
    params: ParamsRequest
    directory: DirectoryObject
    handler: Any


processes: List[Process] = []


def get_process(id) -> Process | None:
    p = [p for p in processes if p if p.id == id]
    if p:
        return p[0]
    return None


class FineTuning:
    def __init__(self):
        self.fs = DirectoryObject(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'finetuning'))

    def get_job(self, id: str) -> Optional[JobResponse]:
        process = get_process(id)
        if not process:
            return None

        return JobResponse(
            id=id,
            name=process.params.name,
            created_at=process.created_at,
            running=process.handler.popen() is None,
            learning_completed=process.learning_completed
        )

    def get_logs(self, id: str) -> Optional[LogsResponse]:
        directory = self.fs.directories.first_or_default(lambda o: o.name == id)
        process = get_process(id)
        if not process:
            return None

        log_output = directory.files.first_or_default(lambda o: o.name == 'output.log')
        errors_output = directory.files.first_or_default(lambda o: o.name == 'error.log')
        docking_output = directory.files.first_or_default(lambda o: o.name == 'docking.log')

        return LogsResponse(
            output=log_output.read_string(),
            docking_output=docking_output.read_string(),
            errors=errors_output.read_string()
        )

    def get_params(self, id: str) -> Optional[ParamsResponse]:
        process = get_process(id)
        if not process:
            return None

        return ParamsResponse(
            center_x=process.params.center_x,
            center_y=process.params.center_y,
            center_z=process.params.center_z,
            size_x=process.params.size_x,
            size_y=process.params.size_y,
            size_z=process.params.size_z,
            batch_size=process.params.batch_size,
            minscore=process.params.minscore,
            epochs=process.params.epochs
        )

    def all_jobs(self) -> List[JobResponse]:
        result = []
        for run_dir in self.fs.children:
            job = self.get_job(run_dir.name)
            if job:
                result.append(job)
        return result

    async def _prepare_pdbqt(self, pdb: UploadFile) -> FileObject:
        tmp = self.fs.add_directory('tmp')
        pdb_file = tmp.add_file(str(uuid.uuid4()) + '.pdb')
        pdb_file.write_bytes(await pdb.read())
        pdbqt = self.fs.add_directory('tmp').add_file(str(uuid.uuid4()) + '.pdbqt')
        prepare_pdbqt = self.fs.files.first_or_default(lambda o: o.name == 'prepare_pdbqt.sh')
        subprocess.run([prepare_pdbqt.full_path, pdb_file.full_path, pdbqt.full_path])
        return pdbqt

    async def prepare_pdbqt(self, pdb: UploadFile) -> FileResponse:
        pdbqt = await self._prepare_pdbqt(pdb)
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
                                 log_file.full_path,
                                 error_file.full_path,
                                 ], shell=False, stdout=sys.stdout, stderr=sys.stderr)

    def prepare_dockstream_config(self,
                                  request: ParamsRequest,
                                  pdbqt: FileObject,
                                  logfile: FileObject,
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
        dockstream_config_json['docking']['header']['logging']['logfile'] = logfile.full_path

        poses = directory.add_file('poses.sdf')
        scores = directory.add_file('scores.sdf')
        dockstream_config_json['docking']['docking_runs'][0]['output']['poses']['poses_path'] = poses.full_path
        dockstream_config_json['docking']['docking_runs'][0]['output']['scores']['scores_path'] = scores.full_path

        print('CONFIG', dockstream_config_json)

        dockstream_config.write_json(dockstream_config_json)

        return dockstream_config

    def prepare_toml(self, minscore: float, batch_size: int, epochs: int, csv_result: FileObject, dockstream_config: FileObject, directory: DirectoryObject) -> FileObject:
        # Configure learning
        reinforcement_learning_config = directory.files.first_or_default(lambda o: o.name == 'RL.toml')
        t = reinforcement_learning_config.read(toml.load, mode='rb')
        t['output_csv'] = directory.add_file('rl_direct.csv').full_path
        t['parameters']['batch_size'] = batch_size
        t['parameters']['summary_csv_prefix'] = csv_result.full_path
        t['diversity_filter']['minscore'] = minscore
        t['stage'][0]['chkpt_file'] = directory.add_file('rl_direct.chkpt').full_path
        t['stage'][0]['max_steps'] = epochs
        t['stage'][0]['scoring']['component'][0]['DockStream']['endpoint'][0] \
            ['params']['configuration_path'] = dockstream_config.full_path

        reinforcement_learning_config.write_string(toml.dumps(t))

        return reinforcement_learning_config

    async def stop(self, id: str):
        process = get_process(id)
        if not process:
            return

        process.handler.kill()

    async def delete(self, id: str):
        global processes
        process = get_process(id)
        if not process:
            return

        process.directory.delete()
        process.handler.kill()
        processes = [p for p in processes if p.id != id]

    async def save_params(self, pdb_file: UploadFile, request: ParamsRequest):
        if not pdb_file:
            raise ReinventException(ErrorCode.PDB_NOT_PROVIDED)

        pdbqt_file = await self._prepare_pdbqt(pdb_file)

        id = str(uuid.uuid4())
        directory = self.fs.copy(self.fs.add_directory(id))
        directory.directories.first_or_default(lambda o: o.name == id).delete()
        pdbqt = directory.add_file(pdbqt_file.name)
        pdbqt.write_bytes(pdbqt_file.read_bytes())

        csv_result = directory.add_file('rl_direct')

        docking_log_file = directory.add_file('docking.log')

        dockstream_config_json = self.prepare_dockstream_config(request, pdbqt, docking_log_file, directory)
        reinforcement_learning_config = self.prepare_toml(request.minscore, request.batch_size, request.epochs, csv_result, dockstream_config_json, directory)

        print('TOML', reinforcement_learning_config.read_string())

        directory.add_file('output.log')
        directory.add_file('error.log')

        processes.append(
            Process(
                id=id,
                running=True,
                name=request.name,
                learning_completed=False,
                created_at=datetime.datetime.utcnow(),
                pdbqt_filename=pdbqt.name,
                pdbqt_content=pdbqt.read_string(),
                params=request,
                directory=directory,
                handler=None
            )
        )

    async def run(self, id: str):
        directory = self.fs.copy(self.fs.add_directory(id))

        reinforcement_learning_config = directory.files.first_or_default(lambda o: o.name == 'RL.toml')
        log_file = directory.add_file('output.log')
        error_file = directory.add_file('error.log')

        process = get_process(id)
        process.handler = self.run_reinforcement_learning(reinforcement_learning_config, log_file, error_file)

    def get_smiles(self, id: str) -> SmilesResponse:
        directory = self.fs.directories.first_or_default(lambda o: o.name == id)
        process = get_process(id)
        if not process:
            return SmilesResponse(smiles=[])

        direct_csv = directory.files.first_or_default(lambda o: o.name == 'rl_direct_1.csv')

        with open(direct_csv.full_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        return SmilesResponse(smiles=[
            Smiles(
                smiles=row['SMILES'],
                drugLikeness=row['QED'],
                score=row['Score']
            ) for row in rows
        ])
