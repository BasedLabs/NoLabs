import os
from urllib.parse import urljoin

import aiohttp
from aiohttp import ClientConnectionError
from flask import Request

from src.server import settings
from src.server.services.experiments_structure_loader import ConformationsExperimentsLoader
from src.server.api_handlers.api_handler import ApiHandler
from src.server.initializers import loggers

import requests


class ConformationsApiClient:
    async def logs(self):
        logs_url = urljoin(settings.CONFORMATIONS_API, 'logs')
        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(logs_url) as resp:
                    return (await resp.json())['get-logs']
        except ClientConnectionError as e:
            loggers.logger.warn('Cannot reach conformations server')
            return

    async def errors(self):
        errors_url = urljoin(settings.CONFORMATIONS_API, 'errors')
        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(errors_url) as resp:
                    return (await resp.json())['conformations-errors']
        except ClientConnectionError as e:
            loggers.logger.warn('Cannot reach conformations server')
            return

    def pipeline(self,
                 pdb_content: str,
                 total_frames: int = 10000,
                 take_frame_every: int = 1000,
                 integrator: str = 'LangevinIntegator',
                 friction_coeff: float = 1,
                 step_size: float = 0.002,
                 temperature_kelvin: float = 273.15,
                 replace_nonstandard_residues: bool = True,
                 add_missing_atoms: bool = True,
                 add_missing_hydrogens: bool = True,
                 ignore_missing: bool = False):
        pipeline_url = urljoin(settings.CONFORMATIONS_API, 'conformations')
        j = {
            'pdb_content': pdb_content,
            'total_frames': total_frames,
            'take_frame_every': take_frame_every,
            'integrator': integrator,
            'friction_coeff': friction_coeff,
            'step_size': step_size,
            'temperature_kelvin': temperature_kelvin,
            'replace_nonstandard_residues': replace_nonstandard_residues,
            'add_missing_atoms': add_missing_atoms,
            'add_missing_hydrogens': add_missing_hydrogens,
            'ignore_missing': ignore_missing
        }
        response = requests.post(pipeline_url, json=j)
        if response.status_code == 200:
            return response.json()['result']

        raise Exception('Cannot obtain simulation result')


class ConformationsApiHandler(ApiHandler):
    def __init__(self):
        self.experiments_loader = ConformationsExperimentsLoader()
        self.api_client = ConformationsApiClient()

    def get_experiments(self):
        return self.experiments_loader.load_experiments()

    def get_experiment(self, request):
        experiment_id = request.args.get('id')
        experiment_name = request.args.get('name')

        if not self.experiments_loader.experiment_exists(experiment_id):
            return {'id': experiment_id, 'name': experiment_name, 'data': {}}

        experiment_data = self.experiments_loader.load_experiment(experiment_id)
        return {'id': experiment_id, 'name': experiment_name, 'data': experiment_data}

    def change_experiment_name(self, request: Request):
        j = request.get_json(force=True)
        experiment_id = j['id']
        experiment_name = j['name']

        self.experiments_loader.rename_experiment(experiment_id, experiment_name)

        return {'status': 200}

    def delete_experiment(self, request: Request):
        j = request.get_json(force=True)
        self.experiments_loader.delete_experiment(j['id'])

        return {'status': 200}

    def inference(self, request: Request) -> dict:
        settings.CONFORMATIONS_IN_PROCESS = True
        try:
            protein_files = request.files.getlist('proteinFileInput')
            experiment_name = request.form['experimentName']
            experiment_id = request.form['experimentId']
            total_frames = int(request.form['totalFramesInput'])
            system_temp = float(request.form['tempInput'])
            take_frame_every = int(request.form['takeFrameEveryInput'])
            step_size = float(request.form['stepSizeInput'])
            replace_nonstandard_residues = 'replaceNonstandardResiduesInput' in request.form
            integrator = request.form['integratorSelect']
            friction_coeff = float(request.form['frictionCoeffInput'])
            add_missing_hydrogens = 'addMissingHydrogensCheckbox' in request.form
            add_missing_atoms = 'addMissingAtomsCheckbox' in request.form
            ignore_missing_atoms = 'ignoreMissingAtomsCheckbox' in request.form

            tmp_pdb = 'tmp.pdb'
            protein_files[0].save(tmp_pdb)
            with open(tmp_pdb, 'r') as f:
                pdb = f.read()
            os.remove(tmp_pdb)

            simulation_result = self.api_client.pipeline(pdb_content=pdb,
                                                         total_frames=total_frames,
                                                         take_frame_every=take_frame_every,
                                                         integrator=integrator,
                                                         friction_coeff=friction_coeff,
                                                         step_size=step_size,
                                                         temperature_kelvin=system_temp,
                                                         replace_nonstandard_residues=replace_nonstandard_residues,
                                                         add_missing_atoms=add_missing_atoms,
                                                         add_missing_hydrogens=add_missing_hydrogens,
                                                         ignore_missing=ignore_missing_atoms)
            if simulation_result:
                self.experiments_loader.store_experiment(experiment_id, simulation_result)
                self.experiments_loader.save_experiment_metadata(experiment_id, experiment_name)

            return {
                'id': experiment_id,
                'name': experiment_name,
                'data': {'pdb': simulation_result}
            }
        finally:
            settings.CONFORMATIONS_IN_PROCESS = False
