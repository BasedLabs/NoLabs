import time
from src.server.services.loaders import DTILoader
from src.server.services.experiment_service import DrugDiscovery, ProteinPropertyPrediction, save_experiment_metadata
from flask import Request
import src.server.services.inference_service as inference_service
from src.server import settings
from src.server.services.fasta_reader import get_sequences
from src.server.services.oboreader import read_obo


class ApiHandler:
    def change_experiment_name(self, request: Request):
        j = request.get_json(force=True)
        experiment_id = j['id']
        experiment_name = j['name']

        return 200


drug_discovery = DrugDiscovery()
protein_prediction = ProteinPropertyPrediction()


class AminoAcidLabApiHandler(ApiHandler):
    def inference(self, request: Request) -> dict:
        amino_acid_input_sequence = request.form['inputSequence']
        amino_acid_input_sequence_files = request.files.getlist('inputSequenceFile')
        experiment_id = request.form['experimentId']
        if not amino_acid_input_sequence:
            amino_acid_input_sequence = \
                [seq for seq in get_sequences(amino_acid_input_sequence_files)][0]
        experiment_id = protein_prediction.run(amino_acid_input_sequence)

        localisation_result = ProteinPropertyPrediction.load_result(experiment_id, 'localisation')
        folding_result = ProteinPropertyPrediction.load_result(experiment_id, 'folding')
        gene_ontology_result = ProteinPropertyPrediction.load_result(experiment_id, 'gene_ontology')
        solubility = inference_service.get_solubility_output(experiment_id, 'solubility')

        obo_graph = read_obo(gene_ontology_result)
        return {'id': experiment_id, 'name': 'PULL IT FROM THE BACK', 'data': {
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

    def get_experiments(self):
        return ProteinPropertyPrediction.load_experiment_names()

    def get_experiment(self, request):
        # get name of the experiment and get EXISTING SAVED data based on this name
        experiment_id = int(request.args.get('id'))

        localisation_result = ProteinPropertyPrediction.load_result(experiment_id, 'localisation')
        folding_result = ProteinPropertyPrediction.load_result(experiment_id, 'folding')
        gene_ontology_result = ProteinPropertyPrediction.load_result(experiment_id, 'gene_ontology')
        solubility = inference_service.get_solubility_output(experiment_id, 'solubility')

        obo_graph = read_obo(gene_ontology_result)
        return {'id': experiment_id, 'name': 'PULL IT FROM THE BACK', 'data': {
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

    def delete_experiment(self):
        # Get name of experiment here and delete it based on this name
        return 200


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

    def delete_experiment(self):
        # Get name of experiment here and delete it based on this name
        return 200


class DrugTargetApiHandler(ApiHandler):
    def inference(self, request):
        ligand_files = request.files.getlist('sdfFileInput')
        protein_files = request.files.getlist('proteinFileInput')
        experiment_name = request.form['experimentName']
        experiment_id = request.form['experimentdId']

        experiment_id = drug_discovery.run(ligand_files=ligand_files, protein_files=protein_files)
        DrugDiscovery.save_experiment_metadata(experiment_id, experiment_name=experiment_name)
        data = DTILoader.get_dti_results(experiment_id)

        return {'id': experiment_id, 'name': experiment_name, 'data': data}

    def get_experiments(self):
        return DrugDiscovery.load_experiment_names()

    def get_experiment(self, request):
        # get name of the experiment and get EXISTING SAVED data based on this name
        experiment_id = request.args.get('id')
        data = DTILoader.get_dti_results(experiment_id)
        return {'id': experiment_id, 'data': data}

    def delete_experiment(self):
        # Get name of experiment here and delete it based on this name
        return 200


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
        time.sleep(10)
        # get name of the experiment and get EXISTING SAVED data based on this name
        experiment_id = int(request.args.get('id'))

        return {'id': experiment_id, 'name': 'Experiment 10', 'data': [{
            'proteinName': "AHAHAHAHAHHAHA2222222222222222",
            'ligandName': 'LALSDLASDLASLDASLDA IAM CRAZYYYY',
            'pdb': open('mock_data/test.pdb').read(),
            'sdf': open('mock_data/test.sdf').read(),
            'affinity': 10
        }]}

    def delete_experiment(self):
        # Get name of experiment here and delete it based on this name
        return 200
