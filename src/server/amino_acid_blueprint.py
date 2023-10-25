from flask import Blueprint
from flask import request

from src.server.api_handlers import AminoAcidLabApiHandler


def resolve_api_endpoints(api_handler: AminoAcidLabApiHandler):
    amino_acid_bp = Blueprint('amino-acid', __name__)

    @amino_acid_bp.route('/inference', methods=['POST'])
    def inference():
        return api_handler.inference(request)

    @amino_acid_bp.route('/experiments')
    def get_experiments():
        return api_handler.get_experiments()

    @amino_acid_bp.route('/load-experiment', methods=['GET'])
    def get_experiment():
        return api_handler.get_experiment(request)

    @amino_acid_bp.route('/delete-experiment', methods=['DELETE'])
    def delete_experiment():
        return api_handler.delete_experiment(request)

    @amino_acid_bp.route('/change-experiment-name', methods=['POST'])
    def change_experiment_name():
        return api_handler.change_experiment_name(request)

    @amino_acid_bp.route('/generate-id', methods=['GET'])
    def generate_id():
        return api_handler.gen_uuid()

    return amino_acid_bp
