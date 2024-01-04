from flask import Blueprint, request

from nolabs.server.api_handlers.amino_acid import ApiHandler

def resolve_api_endpoints(api_handler: ApiHandler):
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

    @amino_acid_bp.route('/load-results', methods=['GET'])
    def get_predictions():
        print("getting predictions")
        return api_handler.get_predictions(request)
    
    @amino_acid_bp.route('/load-experiment-progress', methods=['GET'])
    def get_experiment_progress():
        print("getting experiment progress")
        return api_handler.get_experiment_progress(request)
    
    @amino_acid_bp.route('/load-experiment-instance-progress', methods=['GET'])
    def get_experiment_instance_progress():
        print("getting instance progress")
        return api_handler.get_experiment_instance_progress(request)

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
