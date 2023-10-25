import time
from src.server.services.loaders import DTILoader
from src.server.services.experiment_service import DrugDiscovery, ProteinPropertyPrediction
from flask import Request
import src.server.services.inference_service as inference_service
from src.server import settings
from src.server.services.fasta_reader import get_sequences
from src.server.services.oboreader import read_obo


class ApiHandler:
    pass

drug_discovery = DrugDiscovery(settings.use_gpu, settings.is_test)
protein_prediction = ProteinPropertyPrediction(settings.use_gpu, settings.is_test)


class AminoAcidLabApiHandler(ApiHandler):
    def inference(self, request: Request) -> dict:
        amino_acid_input_sequence = request.form['inputSequence']
        amino_acid_input_sequence_files = request.files.getlist('inputSequenceFile')
        experiment_name = request.form['experimentName']
        experiment_id = request.form['experimentId']
        if not amino_acid_input_sequence:
            amino_acid_input_sequence = \
                [seq for seq in get_sequences(amino_acid_input_sequence_files)][0]
        experiment_id = protein_prediction.run(amino_acid_input_sequence, experiment_id)
        ProteinPropertyPrediction.save_experiment_metadata(experiment_id, experiment_name=experiment_name)
        localisation_result = protein_prediction.load_result(experiment_id, 'localisation')
        folding_result = protein_prediction.load_result(experiment_id, 'folding')
        gene_ontology_result = protein_prediction.load_result(experiment_id, 'gene_ontology')
        solubility_result = protein_prediction.load_result(experiment_id, 'solubility')

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
            'solubility': solubility_result['solubility']
        }}

    def get_experiments(self):
        return ProteinPropertyPrediction.load_experiment_names()

    def get_experiment(self, request):
        # get name of the experiment and get EXISTING SAVED data based on this name
        experiment_id = request.args.get('id')

        if not experiment_id:
            return {}

        localisation_result = protein_prediction.load_result(experiment_id, 'localisation')
        folding_result = protein_prediction.load_result(experiment_id, 'folding')
        gene_ontology_result = protein_prediction.load_result(experiment_id, 'gene_ontology')
        solubility = protein_prediction.load_result(experiment_id, 'solubility')

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

    def change_experiment_name(self, request: Request):
        j = request.get_json(force=True)
        experiment_id = j['id']
        experiment_name = j['name']

        protein_prediction.rename_experiment(experiment_id, experiment_name)
        return {'status': 200}

    def delete_experiment(self, request: Request):
        j = request.get_json(force=True)
        protein_prediction.delete_experiment(j['id'])
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


class DrugTargetApiHandler(ApiHandler):
    def inference(self, request):
        ligand_files = request.files.getlist('sdfFileInput')
        protein_files = request.files.getlist('proteinFileInput')
        experiment_name = request.form['experimentName']
        experiment_id = request.form['experimentId']

        experiment_id = drug_discovery.run(ligand_files=ligand_files, protein_files=protein_files, experiment_id=experiment_id)
        DrugDiscovery.save_experiment_metadata(experiment_id, experiment_name=experiment_name)
        data = drug_discovery.load_result(experiment_id)

        return {'id': experiment_id, 'name': experiment_name, 'data': data}

    def get_experiments(self):
        return DrugDiscovery.load_experiment_names()

    def get_experiment(self, request):
        # get name of the experiment and get EXISTING SAVED data based on this name
        experiment_id = request.args.get('id')

        if not experiment_id:
            return {}

        data = drug_discovery.load_result(experiment_id)
        return {'id': experiment_id, 'data': data}

    def change_experiment_name(self, request: Request):
        j = request.get_json(force=True)
        experiment_id = j['id']
        experiment_name = j['name']

        drug_discovery.rename_experiment(experiment_id, experiment_name)

        return {'status': 200}

    def delete_experiment(self, request: Request):
        j = request.get_json(force=True)
        drug_discovery.delete_experiment(j['id'])

        return {'status': 200}


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
