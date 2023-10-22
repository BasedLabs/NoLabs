from flask import Blueprint
from flask import request

from src.server.api_handlers import DrugTargetApiHandler


def resolve_api_endpoints(api_handler: DrugTargetApiHandler):
    drug_target_bp = Blueprint('drug-target', __name__)

    @drug_target_bp.route('/inference', methods=['POST'])
    def inference():
        return api_handler.inference(request)

    @drug_target_bp.route('/experiments')
    def get_experiments():
        return api_handler.get_experiments()

    @drug_target_bp.route('/load-experiment', methods=['GET'])
    def get_experiment():
        return api_handler.get_experiment(request)

    @drug_target_bp.route('/delete-experiment', methods=['DELETE'])
    def delete_experiment():
        return api_handler.delete_experiment()

    @drug_target_bp.route('/change-experiment-name', methods=['POST'])
    def change_experiment_name():
        return api_handler.change_experiment_name(request)

    return drug_target_bp
