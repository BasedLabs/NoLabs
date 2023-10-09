import argparse

import torch
from flask import Flask, render_template, request
from flask_cors import CORS

from src.server import settings
from src.server.amino_acid_blueprint import amino_acid_bp
from src.server.drug_target_blueprint import drug_target_bp

app = Flask(__name__)
CORS(app)
app.register_blueprint(drug_target_bp, url_prefix='/api/drug-target')
app.register_blueprint(amino_acid_bp, url_prefix='/api/amino-acid')
app.config['TEMPLATES_AUTO_RELOAD'] = True
app.config['FLASK_DEBUG'] = True
parser = argparse.ArgumentParser()
parser.add_argument('--port', default=5000)
parser.add_argument('--host', default='127.0.0.1')
parser.add_argument('--test', default=False)
use_gpu = torch.cuda.is_available()


# TESTS



if __name__ == '__main__':
    args = vars(parser.parse_args())
    is_test = args['test']
    settings.is_test = is_test
    app.run(host=args['host'], port=args['port'])
