import os
from typing import List
from urllib.parse import urljoin

import aiohttp
import requests
from aiohttp import ClientConnectionError
from flask import Request
from werkzeug.datastructures import FileStorage

from src.server import settings
from src.server.api_handlers.api_handler import ApiHandler
from src.server.initializers import loggers
from src.server.services.experiments_structure_loader import ProteinDesignExperimentsLoader


class ProteinDesignApiClient:
    async def logs(self):
        logs_url = urljoin(settings.PROTEIN_DESIGN_IN_PROCESS, 'logs')
        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(logs_url) as resp:
                    return (await resp.json())['get-logs']
        except ClientConnectionError as e:
            loggers.logger.warn('Cannot reach protein design server server')
            return

    async def errors(self):
        errors_url = urljoin(settings.PROTEIN_DESIGN_IN_PROCESS, 'errors')
        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(errors_url) as resp:
                    return (await resp.json())['protein-design-errors']
        except ClientConnectionError as e:
            loggers.logger.warn('Cannot reach protein design server')
            return

    def pipeline(self,
                 pdb_content: str = None,
                 contig: str = '50',
                 symmetry: str = None,
                 timesteps: int = None,
                 hotspots: str = None):
        pipeline_url = urljoin(settings.PROTEIN_DESIGN_API, 'protein-design')
        j = {
            'pdb_content': pdb_content,
            'contig': contig,
            'symmetry': symmetry,
            'timesteps': timesteps,
            'hotspots': hotspots
        }
        response = requests.post(pipeline_url, json=j)
        if response.status_code == 200:
            return response.json()['result']

        raise Exception('Cannot obtain protein design result')


class ProteinDesignApiHandler(ApiHandler):
    def __init__(self):
        self.experiments_loader = ProteinDesignExperimentsLoader()
        self.api_client = ProteinDesignApiClient()

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
        settings.PROTEIN_DESIGN_IN_PROCESS = True
        try:
            pdb_files: List[FileStorage] = request.files.getlist('proteinFileInput')
            contig = request.form['contig']
            symmetry = request.form['symmetry']
            timesteps = int(request.form['timesteps']) if request.form['timesteps'] else None
            hotspots = request.form['hotspots']
            experiment_name = request.form['experimentName']
            experiment_id = request.form['experimentId']

            tmp_pdb = 'tmp.pdb'
            pdb_contents = tuple()
            if pdb_files:
                file_storage = pdb_files[0]
                file_storage.save(tmp_pdb)
                with open(tmp_pdb, 'r') as f:
                    pdb = f.read()
                pdb_contents = (file_storage.name, pdb)
                os.remove(tmp_pdb)

            if not pdb_contents:
                protein_design_result = self.api_client.pipeline(contig=contig,
                                                                 symmetry=symmetry,
                                                                 timesteps=timesteps,
                                                                 hotspots=hotspots)
                self.experiments_loader.store_experiment(experiment_id, protein_design_result)
                self.experiments_loader.save_experiment_metadata(experiment_id, experiment_name)
            else:
                file_name, pdb_content = pdb_contents
                protein_design_result = self.api_client.pipeline(
                    pdb_content=pdb_content,
                    contig=contig,
                    symmetry=symmetry,
                    timesteps=timesteps,
                    hotspots=hotspots)
                self.experiments_loader.store_experiment(experiment_id, protein_design_result, file_name)
                self.experiments_loader.save_experiment_metadata(experiment_id, experiment_name)

            return {
                'id': experiment_id,
                'name': experiment_name,
                'data': {'pdb': protein_design_result}
            }
        finally:
            settings.PROTEIN_DESIGN_IN_PROCESS = False
