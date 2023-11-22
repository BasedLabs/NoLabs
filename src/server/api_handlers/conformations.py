from flask import Request

from src.server.services.conformations.conformations_pipeline import permute_simulation
from src.server.services.experiments_structure_loader import ConformationsExperimentsLoader
from src.server.api_handlers.api_handler import ApiHandler


class ConformationsApiHandler(ApiHandler):
    def __init__(self):
        self.experiments_loader = ConformationsExperimentsLoader()

    def get_experiments(self):
        return self.experiments_loader.load_experiments()

    def get_experiment(self, request):
        experiment_id = request.args.get('id')
        experiment_name = request.args.get('name')

        if not self.experiments_loader.experiment_exists(experiment_id):
            return {'id': experiment_id, 'name': experiment_name, 'data': {}}

        experiment_data = self.experiments_loader.load_experiment(experiment_id)
        return {'id': experiment_id, 'name': experiment_name, 'data':experiment_data}

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
        protein_files = request.files.getlist('proteinFileInput')
        experiment_name = request.form['experimentName']
        experiment_id = request.form['experimentId']

        simulation_result = permute_simulation(protein_files[0])
        if simulation_result:
            print('SIMULATION RESULT', simulation_result)
            self.experiments_loader.store_experiment(experiment_id, simulation_result)
            self.experiments_loader.save_experiment_metadata(experiment_id, experiment_name)

        return {
            'id': experiment_id,
            'name': experiment_name,
            'pdb': simulation_result
        }


