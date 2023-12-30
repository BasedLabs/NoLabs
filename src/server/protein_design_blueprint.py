from flask import Blueprint
from flask import request

from src.server.api_handlers.amino_acid import ApiHandler


def resolve_api_endpoints(api_handler: ApiHandler):
    protein_design_bp = Blueprint('protein-design', __name__)

    @protein_design_bp.route('/inference', methods=['POST'])
    def inference():
        return api_handler.inference(request)

    @protein_design_bp.route('/experiments')
    def get_experiments():
        return api_handler.get_experiments()

    @protein_design_bp.route('/load-experiment', methods=['GET'])
    def get_experiment():
        return api_handler.get_experiment(request)

    @protein_design_bp.route('/delete-experiment', methods=['DELETE'])
    def delete_experiment():
        return api_handler.delete_experiment(request)

    @protein_design_bp.route('/change-experiment-name', methods=['POST'])
    def change_experiment_name():
        return api_handler.change_experiment_name(request)

    @protein_design_bp.route('/generate-id', methods=['GET'])
    def generate_id():
        return api_handler.gen_uuid()

    return protein_design_bp
