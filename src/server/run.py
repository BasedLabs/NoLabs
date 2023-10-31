# import ssl
# ssl._create_default_https_context = ssl._create_unverified_context
import argparse
import torch
import src.server.amino_acid_blueprint as amino_acid_blueprint
import src.server.drug_target_blueprint as drug_target_blueprint

from src.server import settings, factories

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', default=5000)
    parser.add_argument('--host', default='127.0.0.1')
    parser.add_argument('--test', action='store_true')
    args = vars(parser.parse_args())
    is_test = args['test']
    settings.is_test = is_test
    settings.use_gpu = torch.cuda.is_available()
    settings.host = args['host']
    settings.port = args['port']

amino_acid_api_handler, drug_target_api_handler = factories.api_handlers_factory()
from src.server.initializers.initialize_app import init

app, socketio = init()

app.register_blueprint(drug_target_blueprint.resolve_api_endpoints(drug_target_api_handler),
                        url_prefix='/api/drug-target')
app.register_blueprint(amino_acid_blueprint.resolve_api_endpoints(amino_acid_api_handler),
                        url_prefix='/api/amino-acid')
app.config['TEMPLATES_AUTO_RELOAD'] = True
app.config['FLASK_DEBUG'] = True

if __name__ == '__main__':
    print('-- Starting flask server')
    socketio.run(app, host=settings.host, port=settings.port, allow_unsafe_werkzeug=True)
