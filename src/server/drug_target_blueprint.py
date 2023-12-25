from flask import Blueprint
from flask import request
from src.server.api_handlers.drug_target import DrugTargetApiHandler


def resolve_api_endpoints(api_handler: DrugTargetApiHandler):
    drug_target_bp = Blueprint('drug-target', __name__)

    @drug_target_bp.route('/inference', methods=['POST'])
    def inference():
        return api_handler.inference(request)
    
    @drug_target_bp.route('/add-target', methods=['POST'])
    def add_target():
        return api_handler.add_target(request)
    
    @drug_target_bp.route('/load-targets', methods=['GET'])
    def load_targets():
        return api_handler.load_targets(request)
    
    @drug_target_bp.route('/set-binding-pocket', methods=['POST'])
    def store_binding_pocket():
        return api_handler.set_binding_pocket(request)
    
    @drug_target_bp.route('/get-binding-pocket', methods=['GET'])
    def get_binding_pocket():
        return api_handler.get_binding_pocket(request)

    @drug_target_bp.route('/experiments')
    def get_experiments():
        return api_handler.get_experiments()

    @drug_target_bp.route('/load-experiment', methods=['GET'])
    def get_experiment():
        return api_handler.get_experiment(request)
    
    @drug_target_bp.route('/load-results', methods=['GET'])
    def get_predictions():
        print("getting predictions")
        return api_handler.get_predictions(request)
    
    @drug_target_bp.route('/load-experiment-progress', methods=['GET'])
    def get_experiment_progress():
        print("getting experiment progress")
        return api_handler.get_experiment_progress(request)
    
    @drug_target_bp.route('/load-experiment-instance-progress', methods=['GET'])
    def get_experiment_instance_progress():
        print("getting instance progress")
        return api_handler.get_experiment_instance_progress(request)

    @drug_target_bp.route('/delete-experiment', methods=['DELETE'])
    def delete_experiment():
        return api_handler.delete_experiment(request)

    @drug_target_bp.route('/change-experiment-name', methods=['POST'])
    def change_experiment_name():
        return api_handler.change_experiment_name(request)

    @drug_target_bp.route('/generate-id', methods=['GET'])
    def generate_id():
        return api_handler.gen_uuid()

    @drug_target_bp.route('/download-combined-pdb', methods=['POST'])
    def download_combined_pdb():
        return api_handler.download_combined_pdb(request)

    return drug_target_bp
