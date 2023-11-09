from flask import Request

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