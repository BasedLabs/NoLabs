import datetime
import multiprocessing
import os
import re
import subprocess
import time
import tomllib

import toml

import leaf as lf
import uuid
from typing import List, Any, Dict

from fastapi import UploadFile

from leaf import FileObject, DirectoryObject
from leaf.filesystem import ObjectType
from reinvent.exceptions import ReinventException, ErrorCode
from reinvent.api_models import RunFineTuningJobRequest, FineTuningJobResponse, RunInferenceRequest, \
    RunInferenceResponse


class Inference:
    def __init__(self):
        self.finetuning = FineTuning()
        self.fs = lf.DirectoryObject(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'inference'))

    def run(self, request: RunInferenceRequest) -> RunInferenceResponse:
        directory = self.finetuning.fs.directories.first_or_default(lambda o: o.name == request.job_id)
        metadata = directory.files.first_or_default(lambda o: o.name == 'metadata.json').read_json()
        sampling_toml = self.fs.files.first_or_default(lambda o: o.name == 'sampling.toml')
        sampling = sampling_toml.read(tomllib.load)
        output = directory.add_file(f'{request.job_id}.csv')
        sampling['parameters']['chkpt_file'] = directory.files.first_or_default(lambda o: o.name == 'rl_direct.chkpt').full_path
        sampling['parameters']['output_file'] = output.full_path
        new_toml = directory.add_file(f'{request.job_id}.toml')
        new_toml.write_string(toml.dumps(sampling))

        shell = self.fs.files.first_or_default(
            lambda o: o.name == 'sampling.sh')
        subprocess.run([shell.full_path, new_toml.full_path], shell=True, check=True, stdout=subprocess.PIPE, text=True)
        smiles = output.read(tomllib.load)
        return RunInferenceResponse(
            id=request.job_id,
            job_name=metadata['name'],
            smiles=smiles,
            errors=metadata['errors']
        )


class FineTuning:
    def __init__(self):
        self.fs = lf.DirectoryObject(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'finetuning'))

    def get_job(self, id: str) -> FineTuningJobResponse | None:
        metadata_obj = (self.fs
                        .files
                        .first_or_default(lambda o: o.parent.name == id and o.name == 'metadata.json',
                                          recursive=True))

        if not metadata_obj:
            return

        metadata: Dict[str, Any] = metadata_obj.read_json()
        pdb_filename = metadata['pdb_filename']
        pdb_file: lf.FileObject | None = self.fs.files.first_or_default(lambda o: o.name == pdb_filename,
                                                                        recursive=True)
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

    def calculate_progress(self, metadata, current_stage) -> float:
        return float(current_stage) / float(metadata['epochs'])

    def run_refresh_progress_task(self):
        def task():
            while True:
                self.refresh_progress_task()
                time.sleep(60)
        background_process = multiprocessing.Process(target=task)
        background_process.daemon = True
        background_process.start()

    def refresh_progress_task(self):
        directory: DirectoryObject
        for directory in self.fs.where(lambda o: o.type == ObjectType.DIRECTORY):
            metadata_file: FileObject = directory.files.first_or_default(lambda o: o.name == 'metadata.json')
            metadata = metadata_file.read_json()
            log_file: FileObject = directory.files.first_or_default(lambda o: o.name == 'output.log')

            stage = 1.0
            for log in reversed(log_file.read_lines()):
                groups = re.search('Stage [0-9]{1,3}', log)
                if groups:
                    stage = float(re.search('[0-9]{1,3}', groups.group(0)).group(0))
                    break

            errors_file: FileObject = directory.files.first_or_default(lambda o: o.name == 'error.log')
            metadata['errors'] = errors_file.read_string()
            metadata['progress'] = self.calculate_progress(metadata, stage)

            if stage == float(metadata['epochs']):
                metadata['finished_at'] = datetime.datetime.now(datetime.UTC)
                metadata['running'] = False

            metadata_file.write_json(metadata)

    def all_jobs(self) -> List[FineTuningJobResponse]:
        result = []
        for run_dir in self.fs.children:
            job = self.get_job(run_dir.name)
            if job:
                result.append(job)
        return result

    def prepare_pdbqt(self, pdb_file: FileObject, pdbqt: FileObject):
        prepare_pdbqt = self.fs.files.first_or_default(lambda o: o.name == 'prepare_pdbqt.sh')
        subprocess.run([prepare_pdbqt.full_path, pdb_file.full_path, pdbqt.full_path])

    def run_reinforcement_learning(self, RL_toml_file: FileObject, log_file: FileObject, error_file: FileObject) -> int:
        """Returns PID of the process"""
        run_reinforcement_learning_shell = self.fs.files.first_or_default(
            lambda o: o.name == 'start_reinforcement_learning.sh')
        result = subprocess.run([run_reinforcement_learning_shell.full_path,
                                 RL_toml_file.full_path,
                                 log_file.full_path,
                                 error_file.full_path], shell=True, check=True, stdout=subprocess.PIPE, text=True)
        lines = result.stdout.splitlines()
        return int(lines[-1].replace('PID:', '').strip())

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
        t = reinforcement_learning_config.read(tomllib.load)
        t['output_csv'] = directory.add_file('rl_direct.csv').full_path
        t['stage'][0]['chkpt_file'] = directory.add_file('rl_direct.chkpt').full_path
        t['stage'][0]['max_steps'] = request.epochs
        t['stage.scoring.component'][0]['stage.scoring.component.DockStream.endpoint'][0][
            'params.configuration_path'] = dockstream_config.full_path
        reinforcement_learning_config.write_string(toml.dumps(t))

        log_file = directory.add_file('output.log')
        error_file = directory.add_file('error.log')
        pid = self.run_reinforcement_learning(reinforcement_learning_config, log_file, error_file)

        directory.add_file('metadata.json').write_json(
            {
                'name': request.name,
                'running': True,
                'started_at': datetime.datetime.now(datetime.UTC),
                'finished_at': None,
                'progress': 0.0,
                'pdb_filename': pdb_file.filename,
                'errors': '',
                'epochs': request.epochs,
                'pid': pid
            }
        )

        return self.get_job(id)
