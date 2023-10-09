import argparse

import torch
from flask import Flask, render_template, request
from flask_cors import CORS

from src.server import settings
from src.server.amino_acid_blueprint import amino_acid_bp
from src.server.drug_target_blueprint import drug_target_bp

app = Flask(__name__)
CORS(app)
app.register_blueprint(drug_target_bp, url_prefix='/drug-target')
app.register_blueprint(amino_acid_bp, url_prefix='/amino-acid')
app.config['TEMPLATES_AUTO_RELOAD'] = True
app.config['FLASK_DEBUG'] = True
parser = argparse.ArgumentParser()
parser.add_argument('--port', default=5000)
parser.add_argument('--host', default='127.0.0.1')
parser.add_argument('--test', default=False)
use_gpu = torch.cuda.is_available()


@app.route('/')
def index():
    return render_template('index.html')


# TESTS
@app.route('/api/amino-acid-inference', methods=['POST'])
def amino_acid_inference():
    return {
        'sequence': 'AAAAAAAAAA',
        'localisation': {
            'mithochondria': 0.1,
            'nucleus': 0.1,
            'cytoplasm': 0.1,
            'other': 0.1,
            'extracellular': 0.1
        },
        'folding': open('test.pdb').read(),
        'oboGraph': {
            'GO:123': {'name': 'GO:123', 'namespace': 'biological_process', 'edges': {}},
            'GO:234': {'name': 'GO:234', 'namespace': 'biological_process', 'edges': {}}
        },
        'solubility': 0.5
    }


@app.route('/api/drug-target-discovery-inference')
def drug_target_discovery_inference():
    return [{
        'proteinName': "AHAHAHAHAHHAHA",
        'ligandName': 'LALSDLASDLASLDASLDA IAM CRAZYYYY',
        'pdb': open('test.pdb').read(),
        'sdf': open('test.sdf').read(),
        'affinity': 10
    }]


@app.route('/api/drug-target-discovery-experiments')
def drug_target_discovery_experiments():
    return [{
        'id': 'Test1',
        'name': 'Test1'
    }, {
        'id': 'Test2',
        'name': 'Test2'
    }]


@app.route('/api/drug-target-discovery-load-experiment-data/<id>')
def drug_target_discovery_get_experiment():
    experiment_id = request.args.get('id')
    return {
        'name': 'Testtt',
        'id': experiment_id,
        'data': [{
            'proteinName': "AHAHAHAHAHHAHA",
            'ligandName': 'LALSDLASDLASLDASLDA IAM CRAZYYYY',
            'pdb': open('test.pdb').read(),
            'sdf': open('test.sdf').read(),
            'affinity': 10
        }]}


@app.route('/api/drug-target-discovery-add-experiment', methods=['POST'])
def drug_target_discovery_save_experiment():
    return {
        'id': 'Test new from server',
        'name': 'Test new from server'
    }


@app.route('/api/drug-target-discovery-delete-experiment', methods=['DELETE'])
def drug_target_discovery_delete_experiment():
    return 200


if __name__ == '__main__':
    args = vars(parser.parse_args())
    is_test = args['test']
    settings.is_test = is_test
    app.run(host=args['host'], port=args['port'])
