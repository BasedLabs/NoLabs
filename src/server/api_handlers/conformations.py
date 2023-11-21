from flask import Request

from src.server.services.conformations.conformations_pipeline import permute_simulation
from src.server.services.experiments_structure_loader import ConformationsExperimentsLoader
from src.server.api_handlers.api_handler import ApiHandler


class ConformationsApiHandler(ApiHandler):
    def __init__(self):
        self.experiments_loader = ConformationsExperimentsLoader()

    def get_experiments(self):
        d = self.gen_uuid()
        res = {d['id']: 'Test'}
        return res

    def get_experiment(self, request):
        experiment_id = request.args.get('id')
        experiment_name = request.args.get('name')

        return {
            'id': experiment_id,
            'name': experiment_name,
            'data': open('api_handlers/CONFORMATIONS.pdb','r').read()
        }

    def change_experiment_name(self, request: Request):
        return {'result': 'Dont look here! It is a demo'}

    def delete_experiment(self, request: Request):
        return {'result': 'Please no'}

    def inference(self, request: Request) -> dict:
        protein_files = request.files.getlist('proteinFileInput')
        experiment_name = request.form['experimentName']
        experiment_id = request.form['experimentId']

        simulation_result = permute_simulation(pdb_file)
        self.experiments_loader.store_experiment(experiment_id, simulation_result)
        self.experiments_loader.save_experiment_metadata(experiment_id, experiment_name)

        return {
            'id': experiment_id,
            'name': experiment_name,
            'pdb': simulation_result
        }


