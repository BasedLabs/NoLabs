import dataclasses

import datetime
import os
import subprocess
import sys
import csv
from enum import Enum

import toml

import uuid
from typing import List, Any, Optional, Tuple, Set

from fastapi import UploadFile

from leaf import FileObject, DirectoryObject
from starlette.responses import FileResponse

from reinvent.api_models import ConfigurationResponse, ParamsRequest, LogsResponse, ParamsResponse, Smiles, \
    SmilesResponse
from reinvent.exceptions import ReinventException, ErrorCode


class RunType(Enum):
    SAMPLING_SCORING = 1
    RL = 2


@dataclasses.dataclass
class Configuration:
    id: str
    name: str
    created_at: datetime.datetime
    pdbqt_content: str
    pdbqt_filename: str
    params: ParamsRequest
    directory: DirectoryObject
    handler: Any
    learning_started: bool
    generated_smiles_set: Set[str]
    generated_smiles: List[Smiles]

    def __str__(self):
        return f'Process: {self.id}'


configurations: List[Configuration] = []


def get_configuration(config_id) -> Configuration | None:
    p = [p for p in configurations if p if p.id == config_id]
    if p:
        return p[0]
    return None


class Reinvent:
    def __init__(self):
        self.fs = DirectoryObject(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reinvent_configs'))

    def get_config(self, config_id: str) -> Optional[ConfigurationResponse]:
        configuration = get_configuration(config_id)
        if not configuration:
            return None

        running, sampling_allowed = self._check_status(configuration)

        return ConfigurationResponse(
            id=config_id,
            name=configuration.params.name,
            created_at=configuration.created_at,
            running=running,
            sampling_allowed=sampling_allowed
        )

    def _check_status(self, configuration: Configuration) -> Tuple[bool, bool]:
        chkpt = configuration.directory.files.first_or_default(lambda o: o.name == 'rl_direct.chkpt')

        running = configuration.handler is not None and configuration.handler.poll() is None
        sampling_allowed = configuration.learning_started and not running and chkpt is not None

        return (running, sampling_allowed)

    def get_logs(self, config_id: str) -> Optional[LogsResponse]:
        directory = self.fs.directories.first_or_default(lambda o: o.name == config_id)
        config = get_configuration(config_id)
        if not config:
            return None

        log_output = directory.files.first_or_default(lambda o: o.name == 'output.log')
        errors_output = directory.files.first_or_default(lambda o: o.name == 'error.log')
        docking_output = directory.files.first_or_default(lambda o: o.name == 'docking.log')

        return LogsResponse(
            output=log_output.read_string(),
            docking_output=docking_output.read_string(),
            errors=errors_output.read_string()
        )

    def get_params(self, config_id: str) -> Optional[ParamsResponse]:
        config = get_configuration(config_id)
        if not config:
            return None

        return ParamsResponse(
            center_x=config.params.center_x,
            center_y=config.params.center_y,
            center_z=config.params.center_z,
            size_x=config.params.size_x,
            size_y=config.params.size_y,
            size_z=config.params.size_z,
            batch_size=config.params.batch_size,
            minscore=config.params.minscore,
            epochs=config.params.epochs
        )

    def all_configs(self) -> List[ConfigurationResponse]:
        result = []
        for run_dir in self.fs.children:
            config = self.get_config(run_dir.name)
            if config:
                result.append(config)
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

    def run_reinforcement_learning(self, RL_toml_file: FileObject, log_file: FileObject, error_file: FileObject,
                                   directory: DirectoryObject) -> \
            subprocess.Popen[bytes]:
        run_reinforcement_learning_shell = directory.first_or_default(
            lambda o: o.name == 'start_reinforcement_learning.sh')
        return subprocess.Popen([run_reinforcement_learning_shell.full_path,
                                 RL_toml_file.full_path,
                                 error_file.full_path,
                                 log_file.full_path
                                 ], stdout=sys.stdout, stderr=sys.stderr)

    def _prepare_dockstream_config(self,
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

        dockstream_config.write_json(dockstream_config_json)

        return dockstream_config

    def prepare_rl(self, minscore: float, batch_size: int, epochs: int, csv_result: FileObject,
                   dockstream_config: FileObject, directory: DirectoryObject) -> FileObject:
        # Configure learning
        reinforcement_learning_config = directory.files.first_or_default(lambda o: o.name == 'RL.toml')
        t = reinforcement_learning_config.read(toml.load, mode='r')
        t['parameters']['batch_size'] = batch_size
        t['parameters']['summary_csv_prefix'] = csv_result.full_path
        t['diversity_filter']['minscore'] = minscore
        t['stage'][0]['chkpt_file'] = directory.add_file('rl_direct.chkpt').full_path
        t['stage'][0]['max_steps'] = epochs
        t['stage'][0]['scoring']['component'][0]['DockStream']['endpoint'][0] \
            ['params']['configuration_path'] = dockstream_config.full_path

        reinforcement_learning_config.write_string(toml.dumps(t))

        return reinforcement_learning_config

    def prepare_sampling(self, number_of_molecules_to_generate: int, directory: DirectoryObject) -> FileObject:
        chkpt = directory.files.first_or_default(lambda o: o.name == 'rl_direct.chkpt')
        sampling_config: FileObject = directory.files.first_or_default(lambda o: o.name == 'Sampling.toml')
        sampling_output = directory.add_file('sampling_direct.csv')

        sampling_config_toml = sampling_config.read(toml.load, mode='r')
        sampling_config_toml['parameters']['output_file'] = sampling_output.full_path
        sampling_config_toml['parameters']['model_file'] = chkpt.full_path
        sampling_config_toml['parameters']['num_smiles'] = number_of_molecules_to_generate

        sampling_config.write_string(toml.dumps(sampling_config_toml))
        return sampling_config

    def prepare_scoring(self, dockstream_config: FileObject, directory: DirectoryObject) -> FileObject:
        scoring_input = directory.add_file('scoring_input.smi')
        scoring_input.write_string('')

        scoring_config: FileObject = directory.files.first_or_default(lambda o: o.name == 'Scoring.toml')
        scoring_output = directory.add_file('scoring_direct.csv')

        scoring_config_toml = scoring_config.read(toml.load, mode='r')
        scoring_config_toml['output_csv'] = scoring_output.full_path
        scoring_config_toml['parameters']['smiles_file'] = scoring_input.full_path
        scoring_config_toml['scoring']['component'][0]['DockStream']['endpoint'][0] \
            ['params']['configuration_path'] = dockstream_config.full_path

        scoring_config.write_string(toml.dumps(scoring_config_toml))
        return scoring_config

    async def stop(self, config_id: str):
        config = get_configuration(config_id)
        if not config:
            return

        if config.handler:
            config.handler.kill()

    async def delete(self, config_id: str):
        global configurations
        config = get_configuration(config_id)
        if not config:
            return

        config.directory.delete()
        config.handler.kill()
        configurations = [p for p in configurations if p.id != config.id]

    async def save_params(self, pdb_file: UploadFile, request: ParamsRequest):
        if not pdb_file:
            raise ReinventException(ErrorCode.PDB_NOT_PROVIDED)

        await self.delete(request.config_id)

        pdbqt_file = await self._prepare_pdbqt(pdb_file)

        config_id = request.config_id
        directory = self.fs.copy(self.fs.add_directory(config_id))
        directory.directories.first_or_default(lambda o: o.name == config_id).delete()
        pdbqt = directory.add_file(pdbqt_file.name)
        pdbqt.write_bytes(pdbqt_file.read_bytes())

        csv_result = directory.add_file('rl_direct')

        docking_log_file = directory.add_file('docking.log')

        dockstream_config = self._prepare_dockstream_config(request, pdbqt, docking_log_file, directory)
        self.prepare_rl(request.minscore, request.batch_size, request.epochs,
                        csv_result, dockstream_config, directory)
        self.prepare_sampling(10, directory)
        self.prepare_scoring(dockstream_config, directory)

        directory.add_file('output.log')
        directory.add_file('error.log')

        configurations.append(
            Configuration(
                id=config_id,
                name=request.name,
                created_at=datetime.datetime.utcnow(),
                pdbqt_filename=pdbqt.name,
                pdbqt_content=pdbqt.read_string(),
                params=request,
                directory=directory,
                handler=None,
                learning_started=False,
                generated_smiles=[],
                generated_smiles_set=set()
            )
        )

    def run_sampling_and_scoring(self,
                                 number_of_molecules_to_generate: int,
                                 sampling_toml_file: FileObject,
                                 scoring_toml_file: FileObject,
                                 log_file: FileObject, error_file: FileObject, directory: DirectoryObject) -> \
            subprocess.Popen[bytes]:
        """Returns PID of the process"""

        self.prepare_sampling(number_of_molecules_to_generate, directory)

        directory.files.first_or_default(lambda o: o.name == 'sampling_direct.csv').write_string('')
        directory.files.first_or_default(lambda o: o.name == 'scoring_input.smi').write_string('')
        directory.files.first_or_default(lambda o: o.name == 'scoring_direct.csv').write_string('')

        shell = directory.files.first_or_default(
            lambda o: o.name == 'start_sampling.sh')
        return subprocess.Popen([shell.full_path,
                                 sampling_toml_file.full_path,
                                 scoring_toml_file.full_path,
                                 error_file.full_path,
                                 log_file.full_path
                                 ], stdout=sys.stdout, stderr=sys.stderr)

    async def run(self, config_id: str, run_type: RunType, **kwargs):
        config = get_configuration(config_id)

        directory = config.directory

        log_file = directory.files.first_or_default(lambda o: o.name == 'output.log')
        error_file = directory.files.first_or_default(lambda o: o.name == 'error.log')

        if run_type == RunType.RL:
            rl_config = directory.files.first_or_default(lambda o: o.name == 'RL.toml')
            config.handler = self.run_reinforcement_learning(rl_config, log_file, error_file, directory)
        if run_type == RunType.SAMPLING_SCORING:
            sampling_size_key = 'sampling_size'
            if sampling_size_key not in kwargs:
                raise ReinventException(ErrorCode.NUMBER_OF_MOLECULES_REQUIRED)
            sampling_config = directory.files.first_or_default(lambda o: o.name == 'Sampling.toml')
            scoring_config = directory.files.first_or_default(lambda o: o.name == 'Scoring.toml')
            config.handler = self.run_sampling_and_scoring(kwargs[sampling_size_key], sampling_config, scoring_config, log_file, error_file,
                                                           directory)

        config.learning_started = True

    def get_smiles(self, config_id: str) -> SmilesResponse:
        process = get_configuration(config_id)

        if not process:
            return SmilesResponse(
                smiles=[]
            )

        direct = process.directory.files.first_or_default(lambda o: o.name == 'direct.json')

        if not direct:
            direct = process.directory.add_file('direct.json')
        rows = []

        rl_direct = process.directory.files.first_or_default(lambda o: o.name == 'rl_direct_1.csv')

        if rl_direct and rl_direct.size > 0:
            with open(rl_direct.full_path) as f:
                reader = csv.DictReader(f)
                for row in list(reader):
                    smiles = row['SMILES']
                    if smiles not in process.generated_smiles_set:
                        rows.append(
                            {
                                'SMILES': smiles,
                                'drugLikeness': row['QED'],
                                'Score': row['Score'],
                                'Stage': row['step']
                            }
                        )
                        process.generated_smiles.append(
                            Smiles(
                                smiles=smiles,
                                drugLikeness=row['QED'],
                                score=row['Score'],
                                stage=row['step']
                            )
                        )
                        process.generated_smiles_set.add(smiles)

        scoring_direct = process.directory.files.first_or_default(lambda o: o.name == 'scoring_direct.csv')

        if scoring_direct and scoring_direct.size > 0:
            with open(scoring_direct.full_path) as f:
                reader = csv.DictReader(f)
                for row in list(reader):
                    smiles = row['SMILES']
                    if smiles not in process.generated_smiles_set:
                        rows.append(
                            {
                                'SMILES': smiles,
                                'drugLikeness': row['QED'],
                                'Score': row['Score'],
                                'Stage': 'Sampling'
                            }
                        )
                        process.generated_smiles.append(
                            Smiles(
                                smiles=smiles,
                                drugLikeness=row['QED'],
                                score=row['Score'],
                                stage='Sampling'
                            )
                        )
                        process.generated_smiles_set.add(smiles)

        if rows:
            direct.write_json(rows)

        return SmilesResponse(smiles=process.generated_smiles)
