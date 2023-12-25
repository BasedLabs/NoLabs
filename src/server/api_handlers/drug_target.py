import time

from flask import Request, jsonify

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

        return {'id': experiment_id, 'name': experiment_name}

    def add_target(self, request):
        protein_file = request.files['proteinFileInput']
        experiment_name = request.form['experimentName']
        experiment_id = request.form['experimentId']

        self.experiments_loader.save_experiment_metadata(experiment_id, experiment_name=experiment_name)
        self.experiments_loader.store_target(experiment_id, protein_file)

        return {'id': experiment_id, 'name': experiment_name}

    def load_targets(self, request):
        experiment_name = request.args.get('name')
        experiment_id = request.args.get('id')

        targets = self.experiments_loader.load_targets(experiment_id)

        return {'id': experiment_id, 'name': experiment_name, 'targets': targets}

    def set_binding_pocket(self, request):
        experiment_id = request.args.get('id')
        protein_id = request.args.get('proteinId')
        is_pocket_manual = request.args.get('isManualPocket')

        if is_pocket_manual:
            pocket_sequence_ids = request.args.get('pocketIds')
            self.experiments_loader.save_pocket(experiment_id,
                                                protein_id,
                                                is_pocket_manual,
                                                pocket_sequence_ids)
        else:
            self.experiments_loader.save_pocket(experiment_id,
                                                protein_id,
                                                is_pocket_manual)


    def get_binding_pocket(self, request):
        experiment_id = request.args.get('id')
        protein_id = request.args.get('proteinId')
        is_pocket_manual = request.args.get('isManualPocket')

        pocket_ids = []

        pocket_ids = self.experiments_loader.load_pocket(experiment_id,
                                                protein_id,
                                                is_pocket_manual)

        return {"experimentId": experiment_id,
                "proteinId": protein_id,
                "pocketIds": pocket_ids}



    def get_experiments(self):
        return self.experiments_loader.load_experiments()

    def get_experiment(self, request):
        # get name of the experiment and get EXISTING SAVED data based on this name
        experiment_id = request.args.get('id')
        experiment_name = request.args.get('name')

        if not self.experiments_loader.experiment_exists(experiment_id):
            return {'id': experiment_id, 'name': experiment_name, 'data': {}}

        protein_ids = self.experiments_loader.get_protein_ids(experiment_id)

        res = jsonify({'id': experiment_id,
                       'name': experiment_name,
                       'progress': 0,
                       'proteinIds': {protein_id: {'id': protein_id,
                                                'ligandIds': self.experiments_loader.get_ligands_ids(experiment_id=experiment_id, protein_id=protein_id),
                                                  'progress': {'progress': 100.0} } for protein_id in protein_ids} })

        return res

    def get_predictions(self, request):
        # get name of the experiment and get EXISTING SAVED data based on this name
        experiment_id = request.args.get('id')
        experiment_name = request.args.get('name')
        protein_id = request.args.get('proteinId')
        ligand_id = request.args.get('ligandId')

        if not self.experiments_loader.experiment_exists(experiment_id):
            return {'id': experiment_id, 'name': experiment_name, 'data': {}}

        data = self.experiments_loader.load_result(experiment_id, protein_id, ligand_id)
        return {'id': experiment_id, "ligandId": ligand_id, "proteinId": protein_id, 'data': data}

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
        experiment_id = j['experimentId']
        protein_id = j['proteinId']
        ligand_id = j['ligandId']

        data = self.experiments_loader.load_result(experiment_id, protein_id, ligand_id)

        combined_pdb = combine_sdf_pdb(data['sdf'], data['pdb'])

        return {'pdb': combined_pdb}

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