from flask import Request

from server.services.sdf_pdb_combine import combine_sdf_pdb
from src.server.api_handlers.api_handler import ApiHandler
from src.server.services.experiments_structure_loader import DTILabExperimentsLoader


class DrugTargetDemoApiHandler(ApiHandler):
    def __init__(self):
        super().__init__()
        self.experiments_loader = DTILabExperimentsLoader()

    def get_experiments(self):
        return self.experiments_loader.load_experiments()

    def get_experiment(self, request):
        experiment_id = request.args.get('id')
        experiment_name = request.args.get('name')
        data = self.experiments_loader.load_result(experiment_id)
        return {'id': experiment_id, 'data': data}

    def change_experiment_name(self, request: Request):
        return {'result': 'Dont look here! It is a demo'}

    def delete_experiment(self, request: Request):
        return {'result': 'Please no'}

    def inference(self, request: Request) -> dict:
        return {'result': 'This is a demo, deploy docker container or buy us a fancy gpu'}

    def download_combined_pdb(self, request):
        j = request.get_json(force=True)
        experiment_id = j['experiment_id']
        experiment_selected_index = j['selected_index']

        data = self.experiments_loader.load_result(experiment_id)

        combined_pdb = combine_sdf_pdb(data[experiment_selected_index]['sdf'], data[experiment_selected_index]['pdb'])

        return {'pdb': combined_pdb}