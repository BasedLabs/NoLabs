from flask import Request, jsonify

from nolabs.server.services.oboreader import read_obo
from nolabs.server import settings
from nolabs.server.services.experiment_service import ProteinPropertyPrediction
from nolabs.server.api_handlers.api_handler import ApiHandler
from nolabs.server.services.experiments_structure_loader import ProteinLabExperimentsLoader


class AminoAcidLabApiHandler(ApiHandler):
    def __init__(self):

        super().__init__()
        self.experiments_loader = ProteinLabExperimentsLoader()
        self.protein_prediction = ProteinPropertyPrediction(settings.use_gpu, settings.is_test)

    def inference(self, request: Request) -> dict:

        print("FORM: ", request.form)
        amino_acid_input_sequence = request.form['inputSequence']
        fasta_files = request.files.getlist('inputSequenceFile')
        experiment_name = request.form['experimentName']
        experiment_id = request.form['experimentId']

        experiment_id = self.protein_prediction.run(amino_acid_input_sequence, fasta_files, experiment_id)
        self.experiments_loader.save_experiment_metadata(experiment_id, experiment_name=experiment_name)

        return {'id': experiment_id, 'name': experiment_name}

    def get_experiments(self):
        return self.experiments_loader.load_experiments()

    def get_experiment(self, request):
        # get name of the experiment and get EXISTING SAVED data based on this name
        experiment_id = request.args.get('id')
        experiment_name = request.args.get('name')

        if not self.experiments_loader.experiment_exists(experiment_id):
            return {'id': experiment_id, 'name': experiment_name, 'data': {}}

        protein_ids = self.experiments_loader.get_protein_ids(experiment_id)

        print("selecting experiment...")

        res = jsonify({'id': experiment_id, 
                       'name': experiment_name, 
                       'progress': self.experiments_loader.load_experiment_progress(experiment_id),
                       'proteinIds': {protein_id: {'id': protein_id, 'progress': self.experiments_loader.load_protein_progress(experiment_id, protein_id)} for protein_id in protein_ids} })

        return res
    
    def get_experiment_progress(self, request):
        experiment_id = request.args.get('id')
        return {'progress': self.experiments_loader.load_experiment_progress(experiment_id)}
    
    def get_experiment_instance_progress(self, request):
        experiment_id = request.args.get('id')
        protein_id = request.args.get('proteinId')
        print("PROTEIN PROGRESS: ", (experiment_id, protein_id))
        protein_progress = self.experiments_loader.load_protein_progress(experiment_id, protein_id)
        return jsonify({'id': protein_id, 'progress': protein_progress})
    
    def get_predictions(self, request):
        # get name of the experiment and get EXISTING SAVED data based on this name
        experiment_id = request.args.get('id')
        protein_id = request.args.get('proteinId')
        experiment_name = request.args.get('name')

        if not self.experiments_loader.experiment_exists(experiment_id):
            return {'id': experiment_id, 'name': experiment_name, 'data': {}}

        result = self.experiments_loader.load_predictions(experiment_id, protein_id=protein_id)

        localisation_result = result['localisation']
        folding_result = result['folding']
        gene_ontology_result = result['gene_ontology']
        solubility = result['solubility']

        obo_graph = read_obo(gene_ontology_result)
        return {'id': experiment_id, 'name': experiment_name, 'data': {
            'localisation': {
                'mithochondria': localisation_result["Mitochondrial Proteins"],
                'nucleus': localisation_result["Nuclear Proteins"],
                'cytoplasm': localisation_result["Cytosolic Proteins"],
                'other': localisation_result["Other Proteins"],
                'extracellular': localisation_result["Extracellular/Secreted Proteins"]
            },
            'folding': folding_result,
            'oboGraph': obo_graph,
            'solubility': solubility
        }}

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


class AminoAcidLabApiMockHandler(AminoAcidLabApiHandler):
    def inference(self, request: Request) -> dict:
        amino_acid_input_sequence = request.form['inputSequence']
        amino_acid_input_sequence_files = request.files.getlist('inputSequenceFile')
        if request.form['experimentId']:
            experiment_id = int(request.form['experimentId'])
        else:
            experiment_id = None

        return {'id': experiment_id, 'name': 'PULL IT FROM THE BACK', 'data': {
            'sequence': 'AAAAAAAAAA',
            'localisation': {
                'mithochondria': 0.2,
                'nucleus': 0.5,
                'cytoplasm': 0.1,
                'other': 0.4,
                'extracellular': 0.3,
            },
            'folding': open('mock_data/test.pdb').read(),
            'oboGraph': {
                'GO:123': {'name': 'GO:123', 'namespace': 'biological_process', 'edges': {}},
                'GO:234': {'name': 'GO:234', 'namespace': 'biological_process', 'edges': {}}
            },
            'solubility': 0.5
        }}

    def get_experiments(self):
        experiments = [
            {'id': 1, 'name': 'Experiment 13'},
            {'id': 2, 'name': 'Experiment 14'}
        ]
        return experiments

    def get_experiment(self, request):
        # get name of the experiment and get EXISTING SAVED data based on this name
        experiment_id = int(request.args.get('id'))

        return {'id': experiment_id, 'name': 'Test', 'data': {
            'sequence': 'AAAAAAAAAA',
            'localisation': {
                'mithochondria': 0.2,
                'nucleus': 0.5,
                'cytoplasm': 0.1,
                'other': 0.4,
                'extracellular': 0.3,
            },
            'folding': open('mock_data/test.pdb').read(),
            'oboGraph': {
                'GO:123': {'name': 'GO:123', 'namespace': 'biological_process', 'edges': {}},
                'GO:234': {'name': 'GO:234', 'namespace': 'biological_process', 'edges': {}}
            },
            'solubility': 0.5
        }}
