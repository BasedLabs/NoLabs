import time

from flask import Request, jsonify
import numpy as np

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
        experiment_name = request.form['experimentName']
        experiment_id = request.form['experimentId']

        experiment_id = self.drug_discovery.run(experiment_id=experiment_id)
        self.experiments_loader.save_experiment_metadata(experiment_id, experiment_name=experiment_name)

        return {'id': experiment_id, 'name': experiment_name}

    def add_target(self, request):
        protein_file = request.files['proteinFileInput']
        experiment_name = request.form['experimentName']
        experiment_id = request.form['experimentId']

        self.experiments_loader.save_experiment_metadata(experiment_id, experiment_name=experiment_name)
        self.experiments_loader.store_target(experiment_id, protein_file)

        return {'id': experiment_id, 'name': experiment_name}

    def delete_target(self, request):
        data = request.json  # Access data sent in the request body
        experiment_id = data.get('id')
        protein_id = data.get('proteinId')

        self.experiments_loader.delete_target(experiment_id, protein_id)

        return {'id': experiment_id}

    def load_targets(self, request):
        experiment_name = request.args.get('name')
        experiment_id = request.args.get('id')

        targets = self.experiments_loader.load_targets(experiment_id)

        return {'id': experiment_id, 'name': experiment_name, 'targets': targets}

    def predict_3d_structure(self, request):
        experiment_id = request.args.get('id')
        protein_id = request.args.get('proteinId')
        pdb_structure = self.experiments_loader.predict_3d_structure(experiment_id, protein_id)

        return {'pdb': pdb_structure}

    def set_binding_pocket(self, request):
        data = request.json  # Access data sent in the request body
        experiment_id = data.get('id')
        protein_id = data.get('proteinId')
        selected_residues = data.get('selectedResidues', [])
        
        # Convert to NumPy array
        selected_residues_array = np.array(selected_residues)

        self.experiments_loader.set_binding_pocket(experiment_id, 
                                                   protein_id,
                                                   selected_residues_array)

        return {}


    def get_binding_pocket(self, request):
        experiment_id = request.args.get('id')
        protein_id = request.args.get('proteinId')

        pocket_ids = self.experiments_loader.load_binding_pocket(experiment_id,
                                                protein_id)

        return {"experimentId": experiment_id,
                "proteinId": protein_id,
                "pocketIds": pocket_ids}
    
    def predict_binding_pocket(self, request):
        experiment_id = request.args.get('id')
        protein_id = request.args.get('proteinId')

        pocket = self.drug_discovery.predict_pocket(
                                                experiment_id=experiment_id,
                                                protein_id=protein_id
                                                )
        
        return pocket
    
    
    def add_ligand(self, request):
        ligand_file = request.files['ligandFileInput']
        experiment_name = request.form['experimentName']
        experiment_id = request.form['experimentId']

        self.experiments_loader.save_experiment_metadata(experiment_id, experiment_name=experiment_name)
        self.experiments_loader.store_ligand(experiment_id, ligand_file)

        return {'id': experiment_id, 'name': experiment_name}

    def delete_ligand(self, request):
        data = request.json  # Access data sent in the request body
        experiment_id = data.get('id')
        ligand_id = data.get('ligandId')

        self.experiments_loader.delete_ligand(experiment_id, ligand_id)

        return {'id': experiment_id}
    
    def load_ligands(self, request):
        experiment_name = request.args.get('name')
        experiment_id = request.args.get('id')

        ligands = self.experiments_loader.load_ligands(experiment_id)

        return {'id': experiment_id, 'name': experiment_name, 'ligands': ligands}


    def get_experiments(self):
        return self.experiments_loader.load_experiments()

    def get_experiment(self, request):
        experiment_id = request.args.get('id')
        experiment_name = request.args.get('name')

        if not self.experiments_loader.experiment_exists(experiment_id):
            return {}

        return {'id': experiment_id, 'name': experiment_name}
    

    def get_results(self, request):
        experiment_id = request.args.get('id')

        protein_ids = self.experiments_loader.get_protein_ids(experiment_id)

        res = jsonify({'proteinIds': {protein_id: {'id': protein_id,
                                                'ligandIds': self.experiments_loader.get_ligands_ids(experiment_id=experiment_id, protein_id=protein_id),
                                                'ligandResultsAvailable': [self.experiments_loader.check_result_available(experiment_id=experiment_id,
                                                                                                                           protein_id=protein_id, ligand_id=ligand) 
                                                                                                                           for ligand in 
                                                                                                                           self.experiments_loader.get_ligands_ids(experiment_id=experiment_id,
                                                                                                                                                                    protein_id=protein_id)],
                                                'progress': {'progress': 100.0} } for protein_id in protein_ids}})

        return res

    def get_prediction_data(self, request):
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