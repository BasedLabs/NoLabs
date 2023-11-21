from flask import Blueprint
from flask import request

from src.server.api_handlers.api_handler import ApiHandler


def resolve_api_endpoints(api_handler: ApiHandler):
    conformations = Blueprint('conformations', __name__)

    @conformations.route('/inference', methods=['POST'])
    def inference():
        return api_handler.inference(request)

    @conformations.route('/experiments')
    def get_experiments():
        return api_handler.get_experiments()

    @conformations.route('/load-experiment', methods=['GET'])
    def get_experiment():
        return api_handler.get_experiment(request)

    @conformations.route('/delete-experiment', methods=['DELETE'])
    def delete_experiment():
        return api_handler.delete_experiment(request)

    @conformations.route('/change-experiment-name', methods=['POST'])
    def change_experiment_name():
        return api_handler.change_experiment_name(request)

    @conformations.route('/generate-id', methods=['GET'])
    def generate_id():
        return api_handler.gen_uuid()

    return conformations
