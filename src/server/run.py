import argparse

import torch
from flask import Flask, render_template, request
from flask_cors import CORS

from src.server import settings
import src.server.amino_acid_blueprint as amino_acid_blueprint
import src.server.drug_target_blueprint as drug_target_blueprint
from src.server.api_handlers import (AminoAcidLabApiHandler, DrugTargetApiHandler,
                                     AminoAcidLabApiMockHandler, DrugTargetApiMockHandler)

app = Flask(__name__)
CORS(app)
app.register_blueprint(drug_target_blueprint.resolve_api_endpoints(DrugTargetApiMockHandler()),
                       url_prefix='/api/drug-target')
app.register_blueprint(amino_acid_blueprint.resolve_api_endpoints(AminoAcidLabApiMockHandler()),
                       url_prefix='/api/amino-acid')
app.config['TEMPLATES_AUTO_RELOAD'] = True
app.config['FLASK_DEBUG'] = True
parser = argparse.ArgumentParser()
parser.add_argument('--port', default=5000)
parser.add_argument('--host', default='127.0.0.1')
parser.add_argument('--test', default=False)

if __name__ == '__main__':
    args = vars(parser.parse_args())
    is_test = args['test']
    settings.is_test = is_test
    settings.use_gpu = torch.cuda.is_available()
    app.run(host=args['host'], port=args['port'])
