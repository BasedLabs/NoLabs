from flask import Request

from nolabs.server.api_handlers.api_handler import ApiHandler


class ConformationsApiHandler(ApiHandler):
    def __init__(self):
        pass

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
        return {'result': 'This is a demo, deploy docker container or buy us a fancy gpu'}
