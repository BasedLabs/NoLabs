import time

from flask import Request

from src.server.services.sdf_pdb_combine import combine_sdf_pdb
from src.server import settings
from src.server.api_handlers.api_handler import ApiHandler
from src.server.services.experiment_service import DrugDiscovery
from src.server.services.experiments_structure_loader import DTILabExperimentsLoader


class DrugTargetApiHandler(ApiHandler):
    def __init__(self):
        super().__init__()
        self.experiments_loader = DTILabExperimentsLoader()
        self.drug_discovery = DrugDiscovery(settings.use_gpu, settings.is_test)

    def inference(self, request):
        ligand_files = request.files.getlist('sdfFileInput')
        protein_files = request.files.getlist('proteinFileInput')
        experiment_name = request.form['experimentName']
        experiment_id = request.form['experimentId']

        experiment_id = self.drug_discovery.run(ligand_files=ligand_files, protein_files=protein_files,
                                                experiment_id=experiment_id)
        self.experiments_loader.save_experiment_metadata(experiment_id, experiment_name=experiment_name)
        data = self.experiments_loader.load_result(experiment_id)

        return {'id': experiment_id, 'name': experiment_name, 'data': data}

    def get_experiments(self):
        return self.experiments_loader.load_experiments()

    def get_experiment(self, request):
        # get name of the experiment and get EXISTING SAVED data based on this name
        experiment_id = request.args.get('id')
        experiment_name = request.args.get('name')

        if not self.experiments_loader.experiment_exists(experiment_id):
            return {'id': experiment_id, 'name': experiment_name, 'data': {}}

        data = self.experiments_loader.load_result(experiment_id)
        return {'id': experiment_id, 'data': data}

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

    def download_combined_pdb(self, request):
        j = request.get_json(force=True)
        experiment_id = j['experiment_id']
        experiment_selected_index = j['selected_index']

        data = self.experiments_loader.load_result(experiment_id)

        combined_pdb = combine_sdf_pdb(data[experiment_selected_index])

        return {'pdb': 'res'}

class DrugTargetApiMockHandler(DrugTargetApiHandler):
    def inference(self, request):
        ligand_files = request.files.getlist('sdfFileInput')
        protein_files = request.files.getlist('proteinFileInput')
        experiment_id = request.form['experimentId']

        return {'id': experiment_id, 'name': 'PULL IT FROM THE BACK', 'data': [{
            'proteinName': "AHAHAHAHAHHAHA2222222222222222",
            'ligandName': 'LALSDLASDLASLDASLDA IAM CRAZYYYY',
            'pdb': open('mock_data/test.pdb').read(),
            'sdf': open('mock_data/test.sdf').read(),
            'affinity': 10
        }]}

    def get_experiments(self):
        experiments = [
            {'id': 1, 'name': 'Experiment 10'},
            {'id': 2, 'name': 'Experiment 11'}
        ]
        return experiments

    def get_experiment(self, request):
        # get name of the experiment and get EXISTING SAVED data based on this name
        experiment_id = int(request.args.get('id'))

        return {'id': experiment_id, 'name': 'Experiment 10', 'data': [{
            'proteinName': "AHAHAHAHAHHAHA2222222222222222",
            'ligandName': 'LALSDLASDLASLDASLDA IAM CRAZYYYY',
            'pdb': open('mock_data/test.pdb').read(),
            'sdf': open('mock_data/test.sdf').read(),
            'affinity': 10
        }]}

    def download_combined_pdb(self, request):
        pass